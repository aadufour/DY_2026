#! /opt/anaconda3/envs/combine/bin/python3
# Fixed copy of runScans.py
# Change vs original: replaced multiprocessing.Pool with ThreadPoolExecutor
# so that worker threads never try to re-import this script as __main__
# (the root cause of the spawn-bootstrapping crash on macOS / Python 3.12).

from itertools import combinations
import sys
import os
import json
from optparse import OptionParser

def multicore_run(command):
    if getattr(multicore_run, "debug", False):
        print("[DEBUG]", command)
        return
    os.system(command)


point_mode = {
    1: 50,
    2: 1000,
    3: 50
}

secret_options = """--robustFit=1 --setRobustFitTolerance=0.2 --cminDefaultMinimizerStrategy=0 \
--X-rtd=MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=99999999999 --cminFallbackAlgo Minuit2,Migrad,0:0.2 \
--stepSize=0.005 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND """


if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option(
        "--doSplitPoints",
        type="int",
        dest="splitPoints",
        default=0,
        help="How many jobs. Default is 0, which means no splitting",
    )

    parser.add_option(
        "--points",
        type="int",
        dest="points",
        default=None,
        help="User setting of likelihood scan points for GRID algorithm in Combine",
    )

    parser.add_option(
        "--signalPOIs",
        type="str",
        dest="signalPOIs",
        default="",
        help="Comma separated list of parameters to run the scan on --> Define POI of profiled fits",
    )

    parser.add_option(
        "--doOnly",
        type="str",
        dest="doOnly",
        default="",
        help="Comma separated list of pois to consider out of the metadata ones",
    )

    parser.add_option(
        "--verbose",
        type="int",
        dest="verbose",
        default=0,
        help="Verbosity for combine",
    )

    parser.add_option(
        "--stat",
        dest="stat",
        default=False,
        action="store_true",
        help="Freeze all constrained nuisances",
    )

    parser.add_option(
        "--freeze",
        type="str",
        dest="freeze",
        default="r",
        help="Comma separated list of parameters to freeze, by default r",
    )

    parser.add_option(
        "--freezeGroups",
        type="str",
        dest="freezeGroups",
        default="",
        help="Comma separated list of groups and groups to freeze",
    )

    parser.add_option(
        "--unblind",
        dest="unblind",
        default=False,
        action="store_true",
        help="Run unblind fit",
    )

    parser.add_option(
        "--data-asimov",
        dest="dataasim",
        default=False,
        action="store_true",
        help="Run fit with -t -1 --toysFreq",
    )

    parser.add_option(
        "--savePOIs",
        dest="savepois",
        default=False,
        action="store_true",
        help="When running profile fits, save in output the value of the profiled POIs",
    )

    parser.add_option(
        "--track",
        type="str",
        dest="track",
        default="",
        help="Track these parameters during the fit",
    )

    parser.add_option(
        "--randomProf",
        type=int,
        dest="randomprof",
        default=0,
        help="When running profile fits, randomize starting point for profiled pois. Default is 0, no randomizing",
    )

    parser.add_option(
        "--range",
        type="str",
        dest="range",
        default="",
        help="Comma separated list of ranges to be passed to combine",
    )

    parser.add_option(
        "--metadata",
        type="str",
        dest="metadata",
        default="metadata.json",
        help="Metadata file to use e.g. if you want to change ranges. By default metadata.json",
    )

    parser.add_option(
        "--debug",
        dest="debug",
        default=False,
        action="store_true",
        help="Debug mode: print commands only",
    )

    (options, args) = parser.parse_args()

    # mode indicates if we want 1d, 2d or Nd workspaces
    # with N the number of operators
    mode__ = int(sys.argv[1])
    points = point_mode[mode__]
    if options.points != None:
        points = int(options.points)

    if options.splitPoints != 0:
        print(points, options.splitPoints)
        first_point = [i for i in range(0, points, int(points/int(options.splitPoints)))]
        last_points = [i for i in range(int(points/options.splitPoints)-1, points, int(points/int(options.splitPoints)))]


    # initial or scan
    action__ = sys.argv[2]

    with open(options.metadata) as file:
        metadata = json.load(file)

    ops = metadata["operators"].keys()
    if options.doOnly:
        ops = [i for i in options.doOnly.split(",") if i in ops]

    if mode__ == 3: mode__ = len(ops)

    parametrs__ = ""
    if options.signalPOIs:
        parametrs__ = " -P " + " -P ".join([f"k_{op}" for op in options.signalPOIs.split(",")])

    suffix__ = ""
    if len(sys.argv) > 3:
        if sys.argv[3] == "lin":
            suffix__ += "_linear"


    asim = " -t -1 "
    if options.unblind: asim = " "
    if options.dataasim: asim = " -t -1 --toysFreq "

    # track
    track_params = ",".join([i for i in options.track.split(",") if i != ""])
    if track_params != "":
        track_params = " --trackParameters " + track_params + " "

    # freezing
    freeze = ""
    if options.stat: freeze += ",allConstrainedNuisances"
    if options.freeze: freeze += "," + options.freeze
    if options.freezeGroups: freeze += " --freezeNuisanceGroups " + options.freezeGroups
    if freeze.startswith(","): freeze = freeze[1:]
    if freeze != "": freeze = " --freezeParameters " + freeze


    # create a list whose entries will contain the
    # operators for which we want to create the workspace
    combos = list(combinations(ops, mode__))

    cmds = []

    for c in combos:
        print(c)
        operators = ",".join(c)
        name = "_".join(c) + suffix__
        outname = name
        if options.signalPOIs:
            outname = outname + "__" + "_".join(options.signalPOIs.split(","))
        if options.stat:
            outname += "_stat"
        if options.dataasim:
            outname += "_dataasimov"
        if options.unblind:
            outname = outname + "_unblind"
        if mode__ >=3:
            outname += "_profiled"


        setvalue = "--setParameters r=1," + "," + ",".join([f"k_{op}=0" for op in c])
        if options.unblind:
            if not options.freeze:
                setvalue=" "
            elif "r" in options.freeze.split(","):
                setvalue = " --setParameters r=1 "


        redefine = ",".join([f"k_{op}" for op in c])
        if parametrs__ == "":
            pars = " -P " + " -P ".join([f"k_{op}" for op in c])
        else:
            pars = parametrs__

        additional = ""
        if mode__ >= 3:
            if options.savepois:
                additional += f" --saveSpecifiedFunc={redefine} "
            if options.randomprof > 0:
                 additional += f" --pointsRandProf {options.randomprof} "

        print(f"---> Additional: {additional}")
        ranges = ":".join(["k_{}={},{}".format(op, metadata["operators"][op][0], metadata["operators"][op][1]) for op in c])
        if options.range: ranges += ":" + options.range
        if mode__ >= 3:
            ranges = ":".join(["k_{}={},{}".format(op, metadata["operators"][op][0], metadata["operators"][op][1]) for op in c if op in options.signalPOIs.split(",")]) + ":" +  ":".join(["k_{}={},{}".format(op, -300,300) for op in c if op not in options.signalPOIs.split(",")])

        if action__ == "initial":
            add__ = "" if mode__ < 3 else " --floatOtherPOIs=1 "

            if not os.path.isfile(f"model_{name}.root"): continue

            cmd = f"combine -M MultiDimFit model_{name}.root --saveWorkspace -n .initialFit_{outname}  {asim}  --redefineSignalPOIs={redefine} {pars} {setvalue} {freeze} -v {options.verbose} -m 125 {secret_options} {add__} --setParameterRanges={ranges}"
            print(cmd)
            cmds.append(cmd)

        elif action__ == "scan":
            add__ = "" if mode__ < 3 else " --floatOtherPOIs=1 "

            if not os.path.isfile(f"higgsCombine.initialFit_{outname}.MultiDimFit.mH125.root"):
                print("Continue")
                continue

            cmd = f"combineTool.py higgsCombine.initialFit_{outname}.MultiDimFit.mH125.root  -M MultiDimFit --algo grid  -m 125  {asim} --snapshotName MultiDimFit --skipInitialFit --redefineSignalPOIs={redefine} {pars} {setvalue} --setParameterRanges={ranges}  {freeze} -v {options.verbose} --points={points} {secret_options} {add__} {additional} {track_params}"

            if options.splitPoints == 0:
                cmd = cmd + f" -n .{outname}.individual"
                print(cmd)
                cmds.append(cmd)
            else:
                for first,last in zip(first_point, last_points):
                    command = cmd + f" --firstPoint {first} --lastPoint {last} -n .{outname}.individual.POINTS.{first}.{last}"
                    cmds.append(command)


        elif action__ == "singles":
            add__ = "" if mode__ < 3 else " --floatOtherPOIs=1 "

            if not os.path.isfile(f"higgsCombine.initialFit_{outname}.MultiDimFit.mH125.root"):
                print("Continue")
                continue

            cmd = f"combineTool.py higgsCombine.initialFit_{outname}.MultiDimFit.mH125.root  -M MultiDimFit --algo singles  -m 125  {asim} --snapshotName MultiDimFit --redefineSignalPOIs={redefine} {pars} {setvalue} --setParameterRanges={ranges}  {freeze} -v {options.verbose} {secret_options} {add__}"

            cmd = cmd + f" -n .{outname}.individual"
            print(cmd)
            cmds.append(cmd)

    print("Need to launch ", len(cmds), "commands")

    ## Get the number of cores available and run the commands in parallel
    multicore_run.debug = options.debug
    import multiprocessing
    from concurrent.futures import ThreadPoolExecutor  # FIX: threads instead of processes

    cpu_count = multiprocessing.cpu_count()
    num_cores_each_job = min(cpu_count - int(cpu_count/3), len(cmds))

    if options.debug:
        print("DEBUG mode: printing commands, not executing them")
        for cmd in cmds:
            multicore_run(cmd)
    else:
        # FIX: ThreadPoolExecutor uses OS threads, not Python subprocesses.
        # With multiprocessing.Pool (spawn mode, macOS default), each worker
        # re-imports this script as __main__ and hits the Pool() line again,
        # crashing with "bootstrapping phase" RuntimeError.
        # Threads share the process and never trigger that re-import.
        print("Running the processes in multiprocessing mode: {} cores used".format(num_cores_each_job))
        with ThreadPoolExecutor(max_workers=num_cores_each_job) as executor:
            list(executor.map(multicore_run, cmds))

    print("All subprocesses finished")

    if not options.debug and options.splitPoints != 0 and action__ == "scan":
        # Execute the hadd command
        print("Hadding the results!")
        for c in combos:
            name = "_".join(c) + suffix__
            outname = name
            if options.signalPOIs:
                outname = outname + "__" + "_".join(options.signalPOIs.split(","))
            if options.stat:
                outname += "_stat"
            if options.dataasim:
                outname += "_dataasimov"
            if options.unblind:
                outname = outname + "_unblind"
            if mode__ >=3:
                outname += "_profiled"
            command = f"hadd -f higgsCombine.{outname}.individual.MultiDimFit.mH125.root higgsCombine.{outname}.individual.*.MultiDimFit.mH125.root"
            print("Running command: ", command)
            os.system(command)

            # Clean up intermediate root files
            print("Removing intermediate root files...")
            os.system(f"rm -f higgsCombine.{outname}.individual.*.MultiDimFit.mH125.root")
