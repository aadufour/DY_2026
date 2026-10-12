import argparse
import concurrent.futures
import glob
import hashlib
import os
import shutil
from math import ceil
from typing import NewType, Generator

from tqdm import tqdm

from spritz.framework.framework import (  # noqa: F401
    add_dict,
    add_dict_iterable,
    get_analysis_dict,
    get_fw_path,
    read_chunks,
    write_chunks,
    get_batch_cfg,
    low_memory_executor,
)
from spritz.scripts.batch import submit

MERGE_RESULT_FNAME = "tmp_special_"
MERGE_CONDOR_DIR = "merge_condor"

"""
# Result is something like:
{
    "dataset1": {
        # result of single dataset
    }
}
"""
Result = NewType("Result", dict[str, dict])


def read_inputs(inputs: list[str]) -> Generator:
    for input in inputs:
        job_result = read_chunks(input)
        if isinstance(job_result, list):
            for job_result_single in job_result:
                if job_result_single["result"] != {}:
                    yield job_result_single["result"]["real_results"]
        else:
            yield job_result


def check_input(input: Result) -> bool:
    # Returns true if input is ok
    for chunk in input:
        if input["result"] == {} or input["error"] != "":
            return False
    return True


def postprocess_inputs(inputs):
    for input in inputs:
        if MERGE_RESULT_FNAME in input.split("/")[-1]:
            print("removing", input)
            os.remove(input)


def reduction(inputs, reduce_function, output):
    result = reduce_function(read_inputs(inputs))
    postprocess_inputs(inputs)
    write_chunks(result, output)


def split_inputs(inputs, elements_for_task):
    ntasks = ceil(len(inputs) / elements_for_task)
    for i in range(ntasks):
        start = min(i * elements_for_task, len(inputs) - 1)
        stop = min((i + 1) * elements_for_task, len(inputs))
        if start == stop:
            break
        yield slice(start, stop)


def count_tree_tasks(n_inputs, elements_for_task=10):
    """Total number of reduction() calls create_tree() will make across every
    level of the tree, for a progress bar total that covers the whole merge
    (not just one level) up front."""
    total = 0
    n = n_inputs
    while n > elements_for_task:
        n = ceil(n / elements_for_task)
        total += n
    return total + 1  # the final base-case reduction


def create_tree(inputs, reduce_function, output, executor, elements_for_task=10, pbar=None):
    if len(inputs) <= elements_for_task:
        reduction(inputs, reduce_function, output)
        if pbar is not None:
            pbar.update(1)

    else:
        output_dir = "/".join(output.split("/")[:-1])
        output_format = output.split(".")[-1]
        splits = split_inputs(inputs, elements_for_task)
        tasks = []
        new_inputs = []

        futures = {}
        for itask, split in enumerate(split_inputs(inputs, elements_for_task)):
            h = hashlib.new("sha256")
            h.update(str(itask).encode("utf-8"))
            for input in inputs[split]:
                h.update(input.encode("utf-8"))
            h = h.hexdigest()[:10]
            output_tmp = f"{output_dir}/{MERGE_RESULT_FNAME}_{h}.{output_format}"
            future = executor.submit(reduction, inputs[split], reduce_function, output_tmp)
            futures[future] = output_tmp
            new_inputs.append(output_tmp)

        for future in concurrent.futures.as_completed(futures):
            future.result()
            del futures[future]
            if pbar is not None:
                pbar.update(1)

        create_tree(new_inputs, reduce_function, output, executor, elements_for_task, pbar=pbar)


# --split-dataset: 
# every dataset is merged on its own into MERGE_SPLIT_DIR/<dataset>.pkl, so no 
# process ever holds more than one dataset in memory. The result keys of a dataset 
# are never summed with those of another, so this is the same sum as the default
# merge

MERGE_SPLIT_DIR = "results_merged_split"
EFT_KEY_SUFFIXES = ("__eft_combined", "__eft_points", "__eft_covariances")


def split_group(key):
    # EFT datasets store points and covariances under two result keys: keep
    # them together, like the job files do
    for suffix in EFT_KEY_SUFFIXES:
        if key.endswith(suffix):
            return key[: -len(suffix)]
    return key


def _short_hash(*parts):
    h = hashlib.new("sha256")
    for part in parts:
        h.update(str(part).encode("utf-8"))
    return h.hexdigest()[:10]


def pop_results(job_result):
    """Like read_inputs for a single file, but drops every chunk as soon as
    it is over, so a summed chunk is freed right away --> lowers RSS"""
    if isinstance(job_result, list):
        job_result.reverse()
        while job_result:
            chunk = job_result.pop()
            if chunk["result"] != {}:
                yield chunk["result"]["real_results"]
    else:
        yield job_result


def scatter_input(input, tmp_dir):
    """Find the datasets in a job file. A file with a single dataset is used
    as is; a file mixing datasets is rewritten as one piece per dataset.
    Returns {dataset: [file]}"""
    groups = {}
    for real_results in pop_results(read_chunks(input, low_memory=True)):
        for key in list(real_results):
            entry = real_results.pop(key)
            group = groups.setdefault(split_group(key), {})
            group[key] = add_dict(group[key], entry) if key in group else entry

    if len(groups) <= 1:
        return {name: [input] for name in groups}

    pieces = {}
    for name in list(groups):
        piece = f"{tmp_dir}/{MERGE_RESULT_FNAME}_{name}_{_short_hash(input)}.pkl"
        write_chunks(groups.pop(name), piece, low_memory=True)
        pieces[name] = [piece]
    return pieces


def reduction_low_memory(inputs, output):
    """reduction() for the inputs of a single dataset: one input file in
    memory at a time. Returns {result key: sumw}"""
    result = {}
    for input in inputs:
        for real_results in pop_results(read_chunks(input, low_memory=True)):
            result = add_dict(result, real_results)
    write_chunks(result, output, low_memory=True)
    postprocess_inputs(inputs)
    return {key: result[key].get("sumw") for key in result}


def merge_split_dataset(inputs, basepath, cpus, elements_for_task=10):
    split_dir = f"{basepath}/{MERGE_SPLIT_DIR}"
    tmp_dir = f"{split_dir}/tmp"
    old = glob.glob(f"{split_dir}/*.pkl")
    if old:
        print(f"Removing {len(old)} files from a previous split merge in {split_dir}")
        for path in old:
            os.remove(path)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    os.makedirs(tmp_dir)

    inputs_by_group = {}
    sumws = {}
    with low_memory_executor(cpus, preload=["spritz.scripts.merge"]) as executor:
        futures = [executor.submit(scatter_input, input, tmp_dir) for input in inputs]
        for future in tqdm(
            concurrent.futures.as_completed(futures), total=len(futures), desc="indexing"
        ):
            for name, files in future.result().items():
                inputs_by_group.setdefault(name, []).extend(files)
        del futures

        # the same tree as create_tree(), one per dataset, all datasets of a
        # level submitted together so that small datasets fill the workers
        total_tasks = sum(
            count_tree_tasks(len(files), elements_for_task)
            for files in inputs_by_group.values()
        )
        print(f"{len(inputs)} inputs, {len(inputs_by_group)} datasets")
        with tqdm(total=total_tasks, desc="merging") as pbar:
            while inputs_by_group:
                futures = []
                next_level = {}
                for name, files in sorted(inputs_by_group.items()):
                    if len(files) <= elements_for_task:
                        output = f"{split_dir}/{name}.pkl"
                        futures.append(executor.submit(reduction_low_memory, files, output))
                        continue
                    next_level[name] = []
                    for itask, split in enumerate(split_inputs(files, elements_for_task)):
                        output = f"{tmp_dir}/{MERGE_RESULT_FNAME}_{name}_{_short_hash(itask, *files[split])}.pkl"
                        futures.append(executor.submit(reduction_low_memory, files[split], output))
                        next_level[name].append(output)
                for future in concurrent.futures.as_completed(futures):
                    sumws.update(future.result())
                    pbar.update(1)
                inputs_by_group = next_level

    shutil.rmtree(tmp_dir, ignore_errors=True)
    print(sorted(sumws.items()))
    print(f"Merged {len(sumws)} datasets into {split_dir}/")
    if os.path.isfile(f"{basepath}/results_merged_new.pkl"):
        print(
            f"WARNING: {basepath}/results_merged_new.pkl also exists, and spritz-postproc "
            "uses it in preference to the split results: remove it to use these"
        )


def get_args():
    parser = argparse.ArgumentParser(
        description="Merge per-chunk analysis results into results_merged_new.pkl"
    )
    parser.add_argument(
        "--condor",
        action="store_true",
        help="Dispatch the merge to HTCondor instead of running it locally: splits the "
        "raw job_*/chunks_job.pkl into --njobs groups, one condor job per group, each "
        "producing one merged intermediate file under merge_condor/. Run again with "
        "--merge once those finish to produce the final results_merged_new.pkl.",
    )
    parser.add_argument(
        "--merge",
        action="store_true",
        help="Finish a previous --condor dispatch: locally reduce merge_condor/'s "
        "intermediate outputs into results_merged_new.pkl (same local reduction as "
        "the no-flag mode, just starting from far fewer, pre-merged inputs).",
    )
    parser.add_argument(
        "--njobs",
        type=int,
        default=100,
        help="Target number of merged output files when using --condor (default 100) "
        "-- i.e. how many condor jobs to split the raw per-chunk outputs across.",
    )
    parser.add_argument(
        "--split-dataset",
        action="store_true",
        dest="split_dataset",
        help="Merge every dataset on its own into results_merged_split/<dataset>.pkl "
        "instead of one results_merged_new.pkl, to minimize memory: no process holds "
        "more than one dataset and one input file. spritz-postproc reads these when "
        "results_merged_new.pkl is absent. Can be combined with --merge.",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=10,
        help="Number of parallel merge processes (default 10); memory scales with it",
    )
    parser.add_argument(
        "-dr",
        "--dry-run",
        action="store_true",
        dest="dry_run",
        help="With --condor, prepare the job directory but don't actually call "
        "condor_submit (same convention as spritz-batch -dr)",
    )
    return parser.parse_args()


def submit_condor_merge(basepath, njobs, dry_run=False):
    """Dispatch the leaf-level merge work to condor: split the raw
    job_*/chunks_job.pkl files (there can be thousands) into `njobs` groups,
    one condor job per group, each just calling the same reduction()/
    add_dict_iterable() used locally on its own group. This is the
    expensive, embarrassingly-parallel part (per earlier profiling, nearly
    all of a local merge's time is the leaf level); the remaining upper
    tree levels are cheap enough to finish locally afterward with --merge.

    NOTE: workers read their assigned inputs by absolute /gwpool/... path
    directly rather than via HTCondor's transfer_input_files (the shared
    condor_submit() JDL template only supports one fixed file list, not an
    arbitrary per-job set) -- this assumes /gwpool is actually mounted on
    the execute nodes, which is the normal setup for a pool like this one
    but hasn't been verified here. Test with a small --njobs first.
    """
    inputs = [os.path.abspath(p) for p in glob.glob(f"{basepath}/job_*/chunks_job.pkl")]
    if not inputs:
        raise Exception(f"No {basepath}/job_*/chunks_job.pkl files found -- nothing to merge.")

    njobs = min(njobs, len(inputs))
    groups = [inputs[i::njobs] for i in range(njobs)]
    groups = [g for g in groups if len(g) > 0]

    chunks = [
        {
            "data": {"inputs": group, "dataset": f"merge_group_{i}"},
            "weight": 1,
            "result": {},
            "error": "",
        }
        for i, group in enumerate(groups)
    ]

    an_dict = get_analysis_dict()
    batch_config = get_batch_cfg()

    print(f"Merging {len(inputs)} inputs across {len(chunks)} condor jobs -> {len(chunks)} merged files in {MERGE_CONDOR_DIR}/")
    submit(
        chunks,
        path_an=os.path.abspath("."),
        an_dict=an_dict,
        njobs=len(chunks),
        dryRun=dry_run,
        script_name=f"{get_fw_path()}/src/spritz/scripts/merge_worker.py",
        batch_config=batch_config,
        job_dir=MERGE_CONDOR_DIR,
    )
    if dry_run:
        print(f"\nDry run: prepared {MERGE_CONDOR_DIR}/ but did not call condor_submit.")
    else:
        print(
            f"\nSubmitted. Once all jobs in {MERGE_CONDOR_DIR}/ have finished "
            "(check with condor_q / condor_history), run:\n"
            "  spritz-merge --merge\n"
            "to combine them into the final results_merged_new.pkl."
        )


def main():
    args = get_args()
    if args.condor and args.split_dataset:
        raise Exception(
            "--split-dataset applies to the local merge: run 'spritz-merge --condor', "
            "then 'spritz-merge --merge --split-dataset'"
        )
    basepath = os.path.abspath(get_batch_cfg()["BATCH_SYSTEM"])
    output = f"{basepath}/results_merged_new.pkl"
    reduce_function = add_dict_iterable
    elements_for_task = 10
    cpus = args.cpus

    if args.condor:
        submit_condor_merge(basepath, njobs=args.njobs, dry_run=args.dry_run)
        return

    if args.merge:
        inputs = glob.glob(f"{MERGE_CONDOR_DIR}/job_*/chunks_job.pkl")
        if not inputs:
            raise Exception(
                f"No {MERGE_CONDOR_DIR}/job_*/chunks_job.pkl files found -- "
                "did you run 'spritz-merge --condor' first?"
            )
    else:
        inputs = glob.glob(f"{basepath}/job_*/chunks_job.pkl")[:]

    if args.split_dataset:
        inputs = [os.path.abspath(input) for input in inputs]
        merge_split_dataset(inputs, basepath, cpus, elements_for_task=elements_for_task)
        return

    total_tasks = count_tree_tasks(len(inputs), elements_for_task)
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        with tqdm(total=total_tasks, desc="merging") as pbar:
            create_tree(
                inputs,
                reduce_function,
                output,
                executor,
                elements_for_task=elements_for_task,
                pbar=pbar,
            )

    results = read_chunks(output)
    print([(dataset, results[dataset]["sumw"]) for dataset in results])


if __name__ == "__main__":
    main()
