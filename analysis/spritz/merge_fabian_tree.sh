#!/bin/bash
# spritz_fabian tree merge (add_dict_iterable, 10 files per task) with N workers: merge_fabian_tree.sh <config_dir> <N> [seq]
# Run it as a condor job (merge_fabian_tree.jdl: 8-core slot, 4 workers, 120 GB): each worker holds ~20 GB
# with all 235 variations (v3, 2026-11-13); 10 local workers took ~100 GB of the shared login node.
cd "$1"
apptainer exec -B /cvmfs -B /grid_mnt --env SPRITZ_PATH=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian --env PYTHONPATH=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/src:/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif python3 -c "
import sys, concurrent.futures, spritz.scripts.merge as m
cpus = int(sys.argv[1])
import os, glob
from spritz.framework.framework import get_batch_cfg, read_chunks
from spritz.framework.framework import add_dict_iterable
basepath = os.path.abspath(get_batch_cfg()[\"BATCH_SYSTEM\"])
inputs = glob.glob(f\"{basepath}/job_*/chunks_job.pkl\")
output = f\"{basepath}/results_merged_new.pkl\"
from spritz.framework.framework import write_chunks
def reduction(inputs, reduce_function, output):
    # sequential: keep only the running sum + one input in memory (stock reduction loads all inputs at once)
    acc = None
    for inp in inputs:
        for o in m.read_inputs([inp]):
            acc = o if acc is None else reduce_function([acc, o])
    m.postprocess_inputs(inputs)
    write_chunks(acc, output)
if len(sys.argv) > 2 and sys.argv[2] == \"seq\":
    m.reduction = reduction
print(len(inputs), \"inputs ->\", output, \"cpus\", cpus, \"mode\", sys.argv[2:], flush=True)
with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as ex:
    m.create_tree(inputs, add_dict_iterable, output, ex, elements_for_task=10)
r = read_chunks(output)
print(len(r), \"datasets\"); [print(k, r[k][\"sumw\"]) for k in r]
" $2 $3 > merge.log 2>&1
echo MERGE DONE >> merge.log
