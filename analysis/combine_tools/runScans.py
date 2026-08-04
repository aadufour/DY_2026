#!/usr/bin/env python3

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

    parser.add_option(
        "--condor",
        dest="condor",
        default=False,
        action="store_true",
        help="Instead of running commands locally, write per-job LLR T3 HTCondor scripts + a .sub file "
             "under --condor-dir and print the condor_submit command - does NOT submit automatically. "
             "After the condor jobs finish, rerun with --hadd-only to combine the results.",
    )

    parser.add_option(
        "--condor-dir",
        type="str",
        dest="condordir",
        default="condor_jobs",
        help="Directory to write generated condor job scripts / .sub / joblist into (default: condor_jobs)",
    )

    parser.add_option(
        "--cmssw-base",
        type="str",
        dest="cmsswbase",
        default="/grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4",
        help="CMSSW release area to cmsenv into on the worker node (must be visible from LLR worker nodes)",
    )

    parser.add_option(
        "--proxy",
        type="str",
        dest="proxy",
        default="/home/llr/cms/adufour/.t3/proxy.cert",
        help="X509_USER_PROXY path to export in each condor job, mirroring analysis/gridpack conventions",
    )

    parser.add_option(
        "--condor-queue",
        type="str",
        dest="condorqueue",
        default="long",
        help="LLR T3Queue value (default: long, matching analysis/gridpack/launch_lhe_propcorr_all.sub)",
    )

    parser.add_option(
        "--request-memory",
        type="str",
        dest="reqmemory",
        default="2G",
        help="Condor request_memory per job (default: 2G - a combine grid-point fit is much lighter than gridpack generation)",
    )

    parser.add_option(
        "--request-cpus",
        type="int",
        dest="reqcpus",
        default=1,
        help="Condor request_cpus per job (default: 1)",
    )

    parser.add_option(
        "--hadd-only",
        dest="haddonly",
        default=False,
        action="store_true",
        help="Skip building/running scan commands entirely and just hadd whatever "
             "higgsCombine.<name>.individual.POINTS.*.root files are already on disk - "
             "use this after your --condor jobs have finished.",
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

    # if len(sys.argv) > 3:
    #     points = int(sys.argv[3])

    with open(options.metadata) as file:
        metadata = json.load(file)

    ops = list(metadata["operators"].keys())
    if options.doOnly:
        # Filter while preserving metadata.json's own key order, not the order
        # operators were typed in --doOnly - pair names (and therefore
        # model_<name>.root / higgsCombine.<name>...root lookups) are built
        # from combinations() over this list, and the original full batch was
        # built from metadata's order with no --doOnly filtering at all. If
        # this preserved typed order instead, an operator typed early but
        # listed late in metadata.json would get pair names in the opposite
        # order from what already exists on disk, and silently miss them.
        only_set = set(options.doOnly.split(","))
        ops = [i for i in ops if i in only_set]

    if mode__ == 3: mode__ = len(ops)

    parametrs__ = ""
    if options.signalPOIs:
        parametrs__ = " -P " + " -P ".join([f"k_{op}" for op in options.signalPOIs.split(",")])
    #if options.signalPOIs:
    #    ops = options.signalPOIs.split(",")

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

    for c in ([] if options.haddonly else combos):
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

    if options.condor:
        import stat as stat_module

        condordir = os.path.abspath(options.condordir)
        scripts_dir = os.path.join(condordir, "scripts")
        logs_dir = os.path.join(condordir, "logs")
        os.makedirs(scripts_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)
        workdir = os.getcwd()

        joblist_path = os.path.join(condordir, "joblist.txt")
        with open(joblist_path, "w") as jf:
            for i, cmd in enumerate(cmds):
                script_path = os.path.join(scripts_dir, f"job_{i}.sh")
                with open(script_path, "w") as sf:
                    sf.write("#!/bin/bash\n")
                    sf.write("set -e\n")
                    sf.write(f"export X509_USER_PROXY={options.proxy}\n")
                    sf.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
                    sf.write(f"cd {options.cmsswbase}/src\n")
                    sf.write("eval `scramv1 runtime -sh`\n")
                    sf.write(f"cd {workdir}\n")
                    sf.write(cmd + "\n")
                os.chmod(script_path, os.stat(script_path).st_mode | stat_module.S_IEXEC)
                jf.write(script_path + "\n")

        sub_path = os.path.join(condordir, "scan.sub")
        with open(sub_path, "w") as subf:
            subf.write(f"""executable = $(script)
universe = vanilla
output = {logs_dir}/$(Process).out
error = {logs_dir}/$(Process).err
log = {logs_dir}/$(Process).log

request_memory = {options.reqmemory}
request_cpus = {options.reqcpus}

T3Queue = {options.condorqueue}
WNTag = el9
include : /opt/exp_soft/cms/t3/t3queue |
requirements = regexp("llrgrwnvm[0-9]+.in2p3.fr", Machine) == FALSE

priority = -50
max_retries = 1

queue script from {joblist_path}
""")

        print(f"\nWrote {len(cmds)} job scripts to {scripts_dir}/")
        print(f"Wrote joblist to {joblist_path}")
        print(f"Wrote condor submit file to {sub_path}")
        print("\nNOTHING WAS SUBMITTED - review the files above first, especially one of the")
        print("generated job_*.sh scripts (confirm --cmssw-base and --proxy are right for you),")
        print("then test with a small slice (e.g. --doOnly one operator pair) before the full batch.")
        print(f"\nWhen ready:  condor_submit {sub_path}")
        print("Once all jobs finish, rerun this exact command with --hadd-only (instead of --condor) to combine the results.")
        sys.exit(0)

    ## Get the number of cores available and run the commands in parallel
    multicore_run.debug = options.debug

    if options.haddonly:
        print("--hadd-only: skipping local execution, jumping to hadd")
    else:
        import multiprocessing
        from multiprocessing import Pool

        cpu_count = multiprocessing.cpu_count()
        num_cores_each_job = min(cpu_count - int(cpu_count/3), len(cmds))

        if options.debug:
                print("DEBUG mode: printing commands, not executing them")
                for cmd in cmds:
                    multicore_run(cmd)
        else:
            pool = Pool(processes=num_cores_each_job)
            print(("Running the processes in multiprocessing mode: {} cores used".format(num_cores_each_job)))
            pool.map(multicore_run, cmds)
            pool.close()
            pool.join()

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
