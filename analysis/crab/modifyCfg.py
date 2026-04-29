import argparse
import os
import json
import sys
import math
import subprocess

edits = []

Premix ={
    "2023C" :       "/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer23_130X_mcRun3_2023_realistic_v13-v1/PREMIX",
    "2023D" :       "/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer23BPix_130X_mcRun3_2023_realistic_postBPix_v1-v1/PREMIX",
    "2022C" :       "/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX",
    "2022D" :       "/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX",
    "2022E" :       "/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX",
    "2022F" :       "/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX",
    "2022G" :       "/Neutrino_E-10_gun/Run3Summer21PrePremix-Summer22_124X_mcRun3_2022_realistic_v11-v2/PREMIX",
    "2018" :        "/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL18_106X_upgrade2018_realistic_v11_L1v1-v2/PREMIX",
    "2017" :        "/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL17_106X_mc2017_realistic_v6-v3/PREMIX",
    "2016HIPM" :    "/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL16_106X_mcRun2_asymptotic_v13-v1/PREMIX",
    "2016noHIPM":   "/Neutrino_E-10_gun/RunIISummer20ULPrePremix-UL16_106X_mcRun2_asymptotic_v13-v1/PREMIX",
}

class Visitor:
    def __init__(self, modtype):
        self.modtype = modtype
        self.results = []

    def enter(self, obj):
        if isinstance(obj, cms.OutputModule):
            self.results.append(obj)

    def leave(self, obj):
        pass

def UpdateModuleIfExists(mod, attrName, newVal, prefix='process.'):
    if hasattr(mod, attrName):
        attr = getattr(mod, attrName)
        prev = '%s' % attr
        attr.setValue(newVal)
        print('Updating %s from %s to %s' % (attrName, prev, getattr(mod, attrName)))
        edits.append('%s%s.%s = %s' % (prefix, mod, attrName, getattr(mod, attrName)))

def UpdateIfExists(process, modName, attrName, newVal, prefix='process.'):
    if hasattr(process, modName):
        mod = getattr(process, modName)
        if hasattr(mod, attrName):
            attr = getattr(mod, attrName)
            print(attr)
            prev = '%s' % attr
            print(prev)
            newArg = type(attr)(newVal)
            print(newArg)
            attr.setValue(newVal)
            print('Updating %s.%s from %s to %s' % (modName, attrName, prev, getattr(mod, attrName)))
            edits.append('%s%s.%s = %s' % (prefix, modName, attrName, getattr(mod, attrName)))

def Update(process, modName, attrName, newVal, prefix='process.'):
    if hasattr(process, modName):
        mod = getattr(process, modName)
        setattr(mod, attrName, newVal)
        print('Updating %s.%s to %s' % (modName, attrName, getattr(mod, attrName)))
        edits.append('%s%s.%s = %s' % (prefix, modName, attrName, getattr(mod, attrName)))

def SetRandomSeeds(process, newSeed):
    if not hasattr(process, 'RandomNumberGeneratorService'):
        return
    mod = getattr(process, 'RandomNumberGeneratorService')
    for par in mod.parameterNames_():
        UpdateIfExists(mod, par, 'initialSeed', newSeed, prefix='process.RandomNumberGeneratorService.')

def SetEvents(process, nEvents):
    UpdateIfExists(process, 'externalLHEProducer', 'nEvents', nEvents)
    UpdateIfExists(process, 'maxEvents', 'input', nEvents)
    UpdateIfExists(process, 'maxEvents', 'output', nEvents)

def FindOutputModule(process):
    visitor = Visitor(cms.OutputModule)
    for schd in process.schedule:
        schd.visit(visitor)
    if len(visitor.results) == 1:
        return visitor.results[0]
    else:
        raise RuntimeError('Too many OutputModules!')

def SetOutputFileName(process, newName, namedModule=None):
    if namedModule is not None:
        UpdateModuleIfExists(getattr(process, namedModule), 'fileName', 'file:%s' % newName)
    else:
        UpdateModuleIfExists(FindOutputModule(process), 'fileName', 'file:%s' % newName)
def SetInputFileName(process, newName):
    UpdateIfExists(process, 'source', 'fileNames', ['file:%s' % newName])

def SetPremixFiles(process, newList):
    if not hasattr(process, 'mixData'):
        return
    mod = getattr(process, 'mixData')
    if 'input' not in mod.parameterNames_():
        return
    
    UpdateIfExists(mod, 'input', 'fileNames', newList, prefix='process.mixData.')

def checkList(data):
    # Check if dataset is defined in the map. 
    # not the smartest thing to use a list but ok
    if any(data.startswith(year) for year in ["2016", "2017", "2018", "2022", "2023", "2024"]):
        data = Premix[data]

    data = str(data).strip()

    try:
        result = subprocess.check_output(
            ["python3", "get_disk_files.py", "-dataset", data],
            stderr=subprocess.STDOUT 
        )

        AVFile = result.splitlines()
        return AVFile
    except subprocess.CalledProcessError as e:
        print("Error executing: {}".format(e.output))
        return []


def UpdateConfig(inputCfg, outputCfg, events=None, randomSeeds=None, inputFile=None, outputFile=None, gridpack=None, outputModule=None, setLumiOffsets=None):
    # Have to reset sys.argv here in case the inputCfg will do its own VarParsing
    # => this is a bit of a hack!
    sys.argv = [inputCfg]
    exec(compile(open(inputCfg, "rb").read(), inputCfg, 'exec'), globals())
    if events is not None:
        SetEvents(process, int(events))
    if randomSeeds is not None:    
        SetRandomSeeds(process, int(randomSeeds))
    if inputFile is not None:
        SetInputFileName(process, str(inputFile))
    if outputFile is not None:
        SetOutputFileName(process, str(outputFile), namedModule=outputModule)
    if gridpack is not None:
        UpdateIfExists(process, 'externalLHEProducer', 'args', [str(gridpack)], prefix='process.externalLHEProducer.')
    if args.checkPremix is not None:
        availableList = checkList(args.checkPremix)
        SetPremixFiles(process, availableList)

    if setLumiOffsets is not None: 
        eventsPerLumi = int(setLumiOffsets)
        firstLumi = 1 + (int(randomSeeds) - 1) * int(math.ceil(float(events) / float(eventsPerLumi)))
        Update(process, 'source', 'numberEventsInLuminosityBlock', cms.untracked.uint32(eventsPerLumi))
        Update(process, 'source', 'firstLuminosityBlock', cms.untracked.uint32(firstLumi))

    print("strategy is {}".format(args.strategy))
    if args.strategy == 0:
        print("--> Printing to {}".format(outputCfg))
        outFile = open(outputCfg, "w")
        outFile.write(process.dumpPython())
        outFile.close()
    elif args.strategy == 1:
        os.system('cp %s %s' % (inputCfg, outputCfg))
        with open(outputCfg) as ofile:
            cfg = ofile.read()
        cfg = cfg.replace('#{EDITS}', '\n'.join(edits))
        outFile = open(outputCfg, "w")
        outFile.write(cfg)
        outFile.close()



parser = argparse.ArgumentParser()
parser.add_argument('io', nargs=2,
                    help='[input cfg] [output cfg]')
parser.add_argument('--events', type=int, default=None, help='Number of events to process')
parser.add_argument('--randomSeeds', type=int, default=None, help='Set random seeds')
parser.add_argument('--inputFile', type=str, default=None, help='Set input file')
parser.add_argument('--outputFile', type=str, default=None, help='Set output file')
parser.add_argument('--gridpack', type=str, default=None, help='Set gridpack to process')
parser.add_argument('--strategy', type=int, default=0, help='Patching strategy')
parser.add_argument('--setLumiOffsets', type=int, default=None, help='Set this many events per lumiBlock')
parser.add_argument('--outputModule', type=str, default=None, help='Output module to modify')
parser.add_argument('--checkPremix', type=str, default=None, help='If not none, put the dataset or the year you want to check')


args = parser.parse_args()


UpdateConfig(args.io[0], args.io[1], events=args.events, randomSeeds=args.randomSeeds, inputFile=args.inputFile, outputFile=args.outputFile, gridpack=args.gridpack, outputModule=args.outputModule, setLumiOffsets=args.setLumiOffsets)


