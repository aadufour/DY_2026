import argparse
import concurrent.futures
import glob
import hashlib
import os
import sys
from math import ceil
from typing import NewType

from spritz.framework.framework import (  # noqa: F401
    add_dict_iterable,
    read_chunks,
    write_chunks,
    get_batch_cfg
)

parser = argparse.ArgumentParser(description="merge of spritz files")
parser.add_argument("--skip-events", dest="skip_events", action="store_true", default=False,
                    help="Do not merge events fields (saves memory). Only merge histos")
parser.add_argument("--cpus", dest="cpus", type=int, default=16,
                    help="Number of parallel workers (default: 16)")
args, _ = parser.parse_known_args()

MERGE_RESULT_FNAME = "tmp_special_"

"""
# Result is something like:
{
    "dataset1": {
        # result of single dataset
    }
}
"""
Result = NewType("Result", dict[str, dict])
# from typing import TypedDict

# class ChunkResult(TypedDict):


# class Result(TypedDict):
#     results: list[ChunkResult]
#     errors: list[ChunkErred]


def read_inputs(inputs: list[str]) -> list[Result]:
    inputs_obj = []
    for input in inputs:
        job_result = read_chunks(input)
        new_job_result = []
        if isinstance(job_result, list):
            for job_result_single in job_result:
                if job_result_single["result"] != {}:
                    r = job_result_single["result"]["real_results"]
                    if args.skip_events:
                        for ds in r:
                            r[ds].pop("events", None)
                    new_job_result.append(r)
            inputs_obj.extend(new_job_result)
        else:
            inputs_obj.append(job_result)
    return inputs_obj


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
    try:
        inputs_obj = read_inputs(inputs)
        result = reduce_function(inputs_obj)
        postprocess_inputs(inputs)
        print("writing to", output)
        write_chunks(result, output)
    except Exception as e:
        import traceback
        print("ERROR in reduction:", e, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        raise


def split_inputs(inputs, elements_for_task):
    ntasks = ceil(len(inputs) / elements_for_task)
    splits = []
    for i in range(ntasks):
        start = min(i * elements_for_task, len(inputs) - 1)
        stop = min((i + 1) * elements_for_task, len(inputs))
        if start == stop:
            break
        splits.append(slice(start, stop))

    return splits


def create_tree(inputs, reduce_function, output, executor, elements_for_task=10):
    if len(inputs) <= elements_for_task:
        reduction(inputs, reduce_function, output)

    else:
        output_dir = "/".join(output.split("/")[:-1])
        output_format = output.split(".")[-1]
        splits = split_inputs(inputs, elements_for_task)
        tasks = []
        new_inputs = []

        for itask, split in enumerate(splits):
            h = hashlib.new("sha256")
            h.update(str(itask).encode("utf-8"))
            for input in inputs[split]:
                h.update(input.encode("utf-8"))
            h = h.hexdigest()[:10]
            output_tmp = f"{output_dir}/{MERGE_RESULT_FNAME}_{h}.{output_format}"
            tasks.append(
                executor.submit(reduction, inputs[split], reduce_function, output_tmp)
            )
            new_inputs.append(output_tmp)
        concurrent.futures.wait(tasks)
        for task in tasks:
            task.result()

        create_tree(new_inputs, reduce_function, output, executor, elements_for_task)


def main():
    basepath = os.path.abspath(get_batch_cfg()["BATCH_SYSTEM"])
    inputs = glob.glob(f"{basepath}/job_*/chunks_job.pkl")[:]
    output = f"{basepath}/results_merged_new.pkl"
    #print(inputs)
    #print(output)
    reduce_function = sum
    reduce_function = add_dict_iterable
    elements_for_task = 10
    cpus = args.cpus
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        create_tree(
            inputs,
            reduce_function,
            output,
            executor,
            elements_for_task=elements_for_task,
        )

    results = read_chunks(output)
    datasets = results.keys()
    print([(dataset, results[dataset]["sumw"]) for dataset in datasets])


if __name__ == "__main__":
    main()
