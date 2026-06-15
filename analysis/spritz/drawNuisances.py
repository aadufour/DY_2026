#!/usr/bin/env python3

"""
File: drawNuisances.py
Author: Giacomo Boldrini
Date: 2025-02-06
Description: Take a datacard, samples and nuisances, draw the variations only for nuisances defined in the datacard
"""
import matplotlib.pyplot as plt 
import mplhep as hep
hep.style.use("CMS")
from numpy import unique, sqrt

import multiprocessing 
from multiprocessing import Pool


import sys
import json 
from optparse import OptionParser
from tqdm import tqdm
import uproot
import os

from HiggsAnalysis.CombinedLimit.DatacardParser import *
from python.DatacardHelpers import datacardHelper
from python.ShapeFileHelper import shapeFileHelper
from python.DatacardConfigParser import DatacardConfigParser

ratio_fig_style = {
    'figsize': (10, 10),
    'gridspec_kw': {'height_ratios': (3, 1)},
}

shaded_style = {
    'facecolor': (0,0,0,0.3),
    'linewidth': 0
}

NOMINAL_COLOR = "#5790fc"
UP_COLOR = "#f89c20"
DOWN_COLOR = "#e42536"



def plot(d__):
    
    file__ =  d__["file"]
    nominal = d__["nominal"]
    systUp = d__["up"]
    systDown = d__["down"]
    binDC__ = d__["bin"]
    sample = d__["sample"]
    shn = d__["shape"]
    out = d__["out"]
    logy = d__["logy"]

    f = uproot.open(file__) 

    nom = f[nominal]

    up = f[systUp]
    down = f[systDown]

    print(file__, nominal, systUp, systDown)

    nom = nom.to_boost()
    up = up.to_boost()
    down = down.to_boost()

    
    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style)
    fig.subplots_adjust(hspace=.07)  # this controls the margin between the two axes
    

    hep.histplot(nom, ax=ax, color=NOMINAL_COLOR)
    hep.histplot(up, ax=ax, color=UP_COLOR)
    hep.histplot(down, ax=ax, color=DOWN_COLOR)


    ax.set_xlabel("")
    ax.set_ylabel("")
    if logy:
        ax.set_yscale("log")

    
    nom_ratio = nom.copy()
    up_ratio = up.copy()
    down_ratio = down.copy()
    
    nom_ratio.view().value = [nom_ratio.values()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
    
    # nom_ratio.view().variance = [nom_ratio.variances()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
    nom_ratio.view().variance = [0]*len(nom.values())

    up_ratio.view().value = [up_ratio.values()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
    down_ratio.view().value = [down_ratio.values()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1

    # up_ratio.view().variance = [up_ratio.variances()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
    # down_ratio.view().variance = [down_ratio.variances()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1

    up_ratio.view().variance = [0]*len(nom.values()) 
    down_ratio.view().variance = [0]*len(nom.values()) #  if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1

    err_up = [1 + i for i in sqrt(nom_ratio.view().variance)]
    err_down = [1 - i for i in sqrt(nom_ratio.view().variance)]

    err_up.append(err_up[-1])
    err_down.append(err_down[-1])
    

    # rax.fill_between(x=nom_ratio.axes.edges[0], y1=err_down , y2=err_up, step='post', **shaded_style)
    hep.histplot(down_ratio, ax=rax, color=DOWN_COLOR)
    hep.histplot(up_ratio, ax=rax, color=UP_COLOR)
    hep.histplot(nom_ratio, ax=rax, color=NOMINAL_COLOR)

    
    rax.set_ylim(0.7, 1.3)
    rax.set_ylabel('Var / Nom.')
    rax.autoscale(axis='x', tight=True)

    hep.cms.label(loc=0, label="Preliminary", data=True, ax=ax)

    fig.savefig(f"{out}/{binDC__}_{sample}_{shn}.png")
    fig.savefig(f"{out}/{binDC__}_{sample}_{shn}.pdf")
    plt.close(fig)



if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option(
        "--datacard",
        type="string",
        dest="datacard",
        default=None,
        help="The datacard to be read",
    )  

    parser.add_option(
        "--output",
        type="string",
        dest="output",
        default="plots",
        help="Output folder path, by default plots",
    )  

    parser.add_option(
        "--samples",
        type="string",
        dest="samples",
        default=None,
        help="Comma separated list of samples to be plotted, by default everything",
    )  

    parser.add_option(
        "--nuis",
        type="string",
        dest="nuis",
        default=None,
        help="Comma separated list of nuisances to be plotted, by default everything",
    )  

    parser.add_option(
        "--logy",
        action="store_true",
        dest="logy",
        default=False,
        help="Use log scale on the y axis",
    )

    parser.add_option(
        "--ncores",
        type="int",
        dest="ncores",
        default=1,
        help="Number of cores to be used, by default 1",
    )  

    (options, args) = parser.parse_args()   

    print(f"--> Loading datacard {options.datacard}")

    if not os.path.isdir(options.output):
        os.system(f"mkdir {options.output}")
    
    DC = datacardHelper.loadDatacard(options.datacard, parser=parser)
    proc = DC.processes if options.samples == None else options.samples.split(",")
    

    if options.nuis == None:
        nuis = [i[0] for i in DC.systs]
    else:
        #for expression in options.nuis.split(","):
        #    prog = re.compile(expression)
        #    #nuis = [i for i in DC.systs if prog.match(i[0])]
        #    print(nuis[0])
        nuis = options.nuis.split(",")

    fileBinMap = datacardHelper.createfileBinMap(DC, onlyMC=True)

    nuisances_for_samples = {i: datacardHelper.getShapeNuisancesForSample(DC, i) for i in proc}
    for key in nuisances_for_samples.keys():
        print(key)
        nuisances_for_samples[key] = {k:j for k,j in nuisances_for_samples[key].items() if k in nuis}

    print(nuisances_for_samples)

    """
    for file__ in tqdm(fileBinMap.keys()):
        f = uproot.open(file__)
        for sample in proc:
            nuisances_sample = nuisances_for_samples[sample]

            for shn in nuisances_sample.keys():
                for nom_shape, nuis_shape, binDC__ in zip(fileBinMap[file__]["nominalShapes"], fileBinMap[file__]["nuisShapes"],fileBinMap[file__]["datacardBin"]) :
                    if binDC__ not in nuisances_sample[shn]: continue 
                    
                    nominal = nom_shape.replace("$PROCESS", sample).replace("$CHANNEL", binDC__)
                    systUp = nuis_shape.replace(
                        "$PROCESS", sample
                    ).replace("$SYSTEMATIC", shn).replace("$CHANNEL", binDC__) + "Up"

                    systDown = nuis_shape.replace(
                        "$PROCESS", sample
                    ).replace("$SYSTEMATIC", shn).replace("$CHANNEL", binDC__) + "Down"

                    nom = f[nominal]
                    up = f[systUp]
                    down = f[systDown]

                    # print(file__, nominal, systUp, systDown)

                    nom = nom.to_boost()
                    up = up.to_boost()
                    down = down.to_boost()

                    # nom = nom.to_boost()
                    # up = up.to_boost()
                    # down = down.to_boost()
                # 
                    # binedges = nom.axes.edges[0]
                    # bincenters = [(binedges[i]+binedges[i+1])/2 for i in range(len(binedges)-1)]
                # 
                # 
                    # nval = nom.values()
                    # nvar = nom.variances()
                # 
                    # dval = down.values()
                    # dvar = down.variances()
                # 
                    # uval = up.values()
                    # uvar = up.variances()
                    # 
                    # # If i plot more histo with boost it freeezes
                    # # So we convert them to matplotlib histos
                # 
                    # fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style)
                    # fig.subplots_adjust(hspace=.07)  # this controls the margin between the two axes
                # 
                    # #ax.hist(bincenters, bins=binedges, weights=nval,  label=f"Nominal {sample}", color=NOMINAL_COLOR)
                    # #ax.hist(bincenters, bins=binedges, weights=dval,  label=f"{nuis} Down", color=DOWN_COLOR)
                    # #ax.hist(bincenters, bins=binedges, weights=uval,  label=f"{nuis} Up", color=UP_COLOR)
                # 
                    # hep.histplot(nom, ax=ax, label=f"Nominal {sample}", color=NOMINAL_COLOR)
                    # hep.histplot(up, ax=ax, label=f"{nuis} Up", color=UP_COLOR)
                    # print("Done")
                    # hep.histplot(down, ax=ax, label=f"{nuis} Down", color=DOWN_COLOR)
                # 
                    # #ax.legend()
                    # ax.set_xlabel("")
                    # ax.set_ylabel("")
                # 
                    # fig.savefig(f"plots/{binDC__}_{sample}_{shn}.png")
                    # fig.savefig(f"plots/{binDC__}_{sample}_{shn}.pdf")
                    # plt.close(fig)

                    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style)
                    fig.subplots_adjust(hspace=.07)  # this controls the margin between the two axes
                    

                    hep.histplot(nom, ax=ax, label=f"Nominal {sample}", color=NOMINAL_COLOR)
                    hep.histplot(up, ax=ax, label=f"{nuis} Up", color=UP_COLOR)
                    hep.histplot(down, ax=ax, label=f"{nuis} Down", color=DOWN_COLOR)


                    #ax.legend(loc="best")
                    ax.set_xlabel("")
                    ax.set_ylabel("")

                    
                    nom_ratio = nom.copy()
                    up_ratio = up.copy()
                    down_ratio = down.copy()
                    
                    nom_ratio.view().value = [nom_ratio.values()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
                    
                    nom_ratio.view().variance = [nom_ratio.variances()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
                    
                    up_ratio.view().value = [up_ratio.values()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
                    down_ratio.view().value = [down_ratio.values()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1

                    up_ratio.view().variance = [up_ratio.variances()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1
                    down_ratio.view().variance = [down_ratio.variances()[i]/nom.values()[i] if nom.values()[i] != 0 else 1 for i in range(len(nom.values()))] # 1, 1, 1

                    err_up = [1 + i for i in sqrt(nom_ratio.view().variance)]
                    err_down = [1 - i for i in sqrt(nom_ratio.view().variance)]

                    err_up.append(err_up[-1])
                    err_down.append(err_down[-1])
                    

                    rax.fill_between(x=nom_ratio.axes.edges[0], y1=err_down , y2=err_up, step='post', **shaded_style)
                    hep.histplot(down_ratio, ax=rax, color=DOWN_COLOR)
                    hep.histplot(up_ratio, ax=rax, color=UP_COLOR)

                    
                    rax.set_ylim(0.7, 1.3)
                    rax.set_ylabel('Var / Nom.')
                    rax.autoscale(axis='x', tight=True)

                    hep.cms.label(loc=0, label="Preliminary", data=True, ax=ax)

                    fig.savefig(f"plots/{binDC__}_{sample}_{shn}.png")
                    fig.savefig(f"plots/{binDC__}_{sample}_{shn}.pdf")
                    plt.close(fig)
                    
    
    print("Done")

    """
    running = []
    # cycle on the files
    for file__ in fileBinMap.keys():
        f = uproot.open(file__)
        for sample in proc:
            nuisances_sample = nuisances_for_samples[sample]

            for shn in nuisances_sample.keys():
                for nom_shape, nuis_shape, binDC__ in zip(fileBinMap[file__]["nominalShapes"], fileBinMap[file__]["nuisShapes"],fileBinMap[file__]["datacardBin"]) :
                    if binDC__ not in nuisances_sample[shn]: continue 
                    
                    nominal = nom_shape.replace("$PROCESS", sample).replace("$CHANNEL", binDC__)
                    systUp = nuis_shape.replace(
                        "$PROCESS", sample
                    ).replace("$SYSTEMATIC", shn).replace("$CHANNEL", binDC__) + "Up"

                    systDown = nuis_shape.replace(
                        "$PROCESS", sample
                    ).replace("$SYSTEMATIC", shn).replace("$CHANNEL", binDC__) + "Down"
                    
                    d__ = {
                        "file": file__,
                        "bin": binDC__,
                        "sample": sample, 
                        "shape": shn,
                        "nominal": nominal, 
                        "up": systUp,
                        "down": systDown,
                        "out": options.output,
                        "logy": options.logy
                    }

                    running.append(d__)

    pool = Pool(processes=options.ncores) 
    print(("Running the processes in multiprocessing mode: {} cores used".format(options.ncores)))

    pool.map(plot, running)

    # Wait for all processes to finish
    pool.close()
    pool.join()

    print("All subprocesses finished")
    
    

