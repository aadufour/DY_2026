
import concurrent.futures
import json
import os
import subprocess
import sys
from copy import deepcopy

import matplotlib as mpl
import mplhep as hep
import numpy as np
import uproot
import hist
from spritz.framework.framework import get_analysis_dict, get_fw_path

mpl.use("Agg")
from matplotlib import pyplot as plt

d = deepcopy(hep.style.CMS)

d["font.size"] = 12
d["figure.figsize"] = (10, 5)

plt.style.use(d)


def darker_color(color):
    rgb = list(mpl.colors.to_rgba(color)[:-1])
    darker_factor = 4 / 5
    rgb[0] = rgb[0] * darker_factor
    rgb[1] = rgb[1] * darker_factor
    rgb[2] = rgb[2] * darker_factor
    return tuple(rgb)


def make_hist(directory, nuisances, sample):
    
    h = directory[f"histo_{sample}"].to_hist()
    histo = {
        'nom': h.values(),
        'stat_up': h.values() + np.sqrt(h.variances()),
        'stat_down': h.values() - np.sqrt(h.variances()),
    }

    for nuisance in nuisances:
        if nuisances[nuisance]["type"] == "rateParam":
            continue
        if nuisances[nuisance]["type"] == "auto":
            continue
        if nuisances[nuisance]["type"] == "stat":
            continue
        name = nuisances[nuisance]["name"]

        if sample not in nuisances[nuisance]["samples"]:
            histo[f"{nuisance}_up"] = h.values().copy()
            histo[f"{nuisance}_down"] = h.values().copy()
            continue

        if nuisances[nuisance]["type"] == "lnN":
            scaling = float(nuisances[nuisance]["samples"][sample])
            histo[f"{nuisance}_up"] = scaling * h.values()
            histo[f"{nuisance}_down"] = 1.0 / scaling * h.values()
        else:
            histo[f"{nuisance}_up"] = directory[f"histo_{sample}_{name}Up"].values().copy()
            histo[f"{nuisance}_down"] = directory[f"histo_{sample}_{name}Down"].values().copy()
    
    return histo


def vary_systematics(histos, samples, nuisances):
    v_syst_bkg = {
        syst: {
            "up": np.zeros_like(histos[samples[0]]['nom']), 
            "down": np.zeros_like(histos[samples[0]]['nom']),
        }
        for syst in nuisances
    }

    for sample in samples:
        for nuisance in nuisances:
            if sample not in nuisances[nuisance]["samples"]:
                v_syst_bkg[nuisance]["up"] += histos[sample]["nom"].copy()
                v_syst_bkg[nuisance]["down"] += histos[sample]["nom"].copy()
            else:
                v_syst_bkg[nuisance]["up"] += histos[sample][nuisance + "_up"].copy()
                v_syst_bkg[nuisance]["down"] += histos[sample][nuisance + "_down"].copy()

    return v_syst_bkg


def plot_data(ax, x, y, zorder, color='black', label='Data', divide=None):
    if divide is None:
        divide = np.ones_like(x)
    ax.errorbar(
        x,
        y['nom']/divide,
        yerr=(y['down']/divide, y['up']/divide),
        fmt="o",
        markersize=4,
        label=f"{label} [{int(round(np.sum(y['nom']), 0))}]",
        zorder=zorder,
        color=color
    )

def plot_mc_stack(ax, edges, histos, colors, divide=None):
    if divide is None:
        divide = np.ones_like(edges[1:])
    for i,sample in enumerate(histos.keys()):
        y = histos[sample]['nom'].copy()
        ax.stairs(
            (y+sum([histos[sample]['nom'] for sample in list(histos.keys())[:i]]))/divide,
            edges,
            label=sample + f" [{int(round(np.sum(y), 0))}]",
            fill=True,
            zorder=-i,
            color=colors[sample],
            edgecolor=darker_color(colors[sample]),
            linewidth=1.0,
        )

def plot_mc_tot(ax, edges, y, zorder, label='Tot MC', unc=True, fill=False, color='darkgrey', divide=None):
    if divide is None:
        divide = np.ones_like(edges[1:])
    if unc:
        unc_up = round(np.sum(y['up']) / np.sum(y['nom']) * 100, 2)
        unc_down = round(np.sum(y['down']) / np.sum(y['nom']) * 100, 2)
        ax.stairs(
            (y['nom']+y['up'])/divide,
            edges,
            baseline=(y['nom']-y['down'])/divide,
            label=f"Syst [-{unc_down}, +{unc_up}]%",
            fill=True,
            hatch="///",
            color='black',
            facecolor="none",
            zorder=zorder+1,
        )

    ax.stairs(
        y['nom']/divide, 
        edges, 
        label=f"{label} [{int(round(np.sum(y['nom']), 0))}]", 
        color=color,
        fill=fill,
        zorder=zorder
    )

def plot_ratio(ax, x, edges, numerator, denominator, yrange=None):
    denominator['nom'] = np.where(denominator['nom'] >= 1e-6, denominator['nom'], 1e-6)
    # plot denominator
    ax.stairs(
        1+denominator['up']/denominator['nom'],
        edges,
        baseline=1-denominator['down']/denominator['nom'],
        fill=True,
        color=denominator.get('color','black'),
        alpha=0.25,
        zorder=-2*(len(numerator)+1)
    )
    ax.plot(
        edges, 
        np.ones_like(edges), 
        color=denominator.get('color','black'), 
        linestyle="dashed",
        zorder=-2*len(numerator)-1
    )

    # plot ratio
    if len(numerator)>1:
        offsets = np.linspace(-0.3, 0.3, len(numerator))
    else:
        offsets = np.array([0.0])
    for i in range(len(numerator)):
        edge, x_i = 0., x.copy()
        for j in range(len(x)):
            width = x[j]-edge
            edge = x[j]+width
            x_i[j] = x[j]+offsets[i]*width
        if numerator[i].get('is_data'):
            ax.errorbar(
                np.array(x_i),
                numerator[i]['nom']/denominator['nom'],
                yerr=(numerator[i]['down']/denominator['nom'], numerator[i]['up']/denominator['nom']),
                fmt="o",
                markersize=4,
                color=numerator[i]['color'],
                zorder=-2*len(numerator)
            )
        else:
            ax.stairs(
                (numerator[i]['nom']+numerator[i]['up'])/denominator['nom'],
                edges,
                baseline=(numerator[i]['nom']-numerator[i]['down'])/denominator['nom'],
                fill=True, 
                color=numerator[i]['color'],
                alpha=0.25,
                zorder=-2*len(numerator)+i,
            )
            ax.stairs(
                numerator[i]['nom']/denominator['nom'], 
                edges, 
                #label=f"{label} [{int(round(np.sum(y), 0))}]", 
                #color=numerator[i]['color'],
                fill=False,
                edgecolor=numerator[i]['color'],
                linewidth=1.0,
                zorder=-2*len(numerator)+1+i,
            )
    
    ymin = min([
        np.min(
            np.where(
                numerator[i]['nom']/denominator['nom']>0,
                numerator[i]['nom']/denominator['nom'],
                np.ones_like(numerator[i]['nom'])
            )
        ) for i in range(len(numerator))
    ]) if yrange is None else yrange[0]
    ymax = max([
        np.max(numerator[i]['nom']/denominator['nom'])
        for i in range(len(numerator))
    ]) if yrange is None else yrange[1]
    
    ylim = min(1., 1.15*max(abs(1-ymax), abs(1-ymin)))
    #ax.set_ylim(min(1-ylim, 0.94), max(1+ylim, 1.06))
    ax.set_xlim(np.min(edges), np.max(edges))

    ax.set_ylim(0.91, 1.09)


def plot(
    input_file,
    region,
    variable,
    weights,
    samples,
    nuisances,
    lumi,
    colors,
    year_label,
    variable_label=None,
):
    print("Doing ", region, variable)

    directory = input_file[f"{region}/{variable}/nominal"]
    mc_samples = [x for x in samples if not samples[x].get('is_data', False)]
    bkg_samples = [x for x in mc_samples if not samples[x].get('is_signal', False)]
    
    # get the histograms
    histos = {
        sample: make_hist(directory, nuisances, sample)
        for sample in samples
    }

    # prepare total MC histogram
    v_syst = vary_systematics(histos, mc_samples, nuisances)
    ymc = sum([histos[x]['nom'] for x in mc_samples])
    ymc_up = np.sqrt(sum([np.square(v_syst[x]['up'].copy() - ymc) for x in v_syst]))
    ymc_down = np.sqrt(sum([np.square(v_syst[x]['down'].copy() - ymc) for x in v_syst]))

    # blinding
    #signal_tot = ymc - ymc_bkg
    #significance = signal_tot / ymc_bkg
    #blind_mask = significance > 0.10
    blind_mask = np.zeros_like(ymc) # not used right now

    # prepare data histogram
    if 'Data' in histos:
        ydata = histos['Data']['nom'].copy()
        ydata_up = histos['Data']['stat_up'].copy() - ydata
        ydata_down = ydata - histos['Data']['stat_down'].copy()
    else:
        ydata = np.zeros_like(ymc)
        ydata_up = np.zeros_like(ymc)
        ydata_down = np.zeros_like(ymc)
    if 'sr' in region:
        ydata = np.where(blind_mask, np.zeros_like(ydata), ydata)
        ydata_up = np.where(blind_mask, np.zeros_like(ydata), ydata_up)
        ydata_down = np.where(blind_mask, np.zeros_like(ydata), ydata_down)

    # set up figure
    axis = directory[f"histo_{mc_samples[0]}"].to_hist().axes[0]
    x = axis.centers
    edges = axis.edges
    widths = axis.widths
    variable_binwidth = isinstance(axis, hist.axis.Variable)

    fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3,1]}, dpi=200)
    fig.tight_layout(pad=-0.5)
    hep.cms.label(region, data=True, lumi=round(lumi, 2), ax=ax[0], year=year_label)

    # plot data
    plot_data(
        ax[0], 
        x, 
        {'nom': ydata, 'up': ydata_up, 'down': ydata_down}, 
        zorder=2,
        divide=widths if variable_binwidth else None
    )

    # plot MC
    plot_mc_stack(
        ax[0], 
        edges, 
        {k:v for k,v in histos.items() if k in mc_samples}, 
        colors,
        divide=widths if variable_binwidth else None
    )

    # plot total MC and uncertainty
    plot_mc_tot(
        ax[0], 
        edges, 
        {'nom': ymc, 'up': ymc_up, 'down': ymc_down}, 
        zorder=1,
        divide=widths if variable_binwidth else None
    )

    # finalize upper panel
    ax[0].set_yscale("log")
    ax[0].legend(
        loc="upper center",
        frameon=True,
        ncols=3,
        framealpha=0.8,
        fontsize=8,
    )
    ax[0].set_ylim(
        max(1e-4 if variable_binwidth else 0.5, np.min(0.1*ydata/widths if variable_binwidth else histos[mc_samples[0]]['nom'])), 
        max(np.max(ydata/widths if variable_binwidth else ydata), np.max(ymc/widths if variable_binwidth else ymc)) * 5e2
    )
    ax[0].set_ylabel("Events/(bin width)" if variable_binwidth else "Events")

    # lower panel
    plot_ratio(
        ax[1], 
        x, 
        edges, 
        [{'nom': ydata, 'up': ydata_up, 'down': ydata_down, 'is_data': True, 'color': 'black'}], 
        {'nom': ymc, 'up': ymc_up, 'down': ymc_down}
    )

    ax[1].set_ylabel("DATA / MC")
    if variable_label:
        ax[1].set_xlabel(variable_label)
    else:
        ax[1].set_xlabel(variable)

    fig.savefig(
        f"plots/{region}_{variable}.png",
        facecolor="white",
        pad_inches=0.1,
        bbox_inches="tight",
    )

    plt.close()

    for weight in weights:
        directory = input_file[f"{region}/{variable}/{weight}"]
        histos_weight = {
            sample: make_hist(directory, nuisances, sample)
            for sample in mc_samples
        }

        # prepare total MC histogram
        v_syst = vary_systematics(histos_weight, mc_samples, nuisances)
        ymc_weight = sum([histos_weight[x]['nom'] for x in mc_samples])
        ymc_up_weight = np.sqrt(sum([np.square(v_syst[x]['up'].copy() - ymc_weight) for x in v_syst]))
        ymc_down_weight = np.sqrt(sum([np.square(v_syst[x]['down'].copy() - ymc_weight) for x in v_syst]))

        # set up figure
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3,2]}, dpi=200)
        fig.tight_layout(pad=-0.5)
        hep.cms.label(region, data=True, lumi=round(lumi, 2), ax=ax[0], year=year_label)

        # plot
        plot_data(
            ax[0], 
            x,
            {'nom': ydata, 'up': ydata_up, 'down': ydata_down}, 
            zorder=0,
            divide=widths if variable_binwidth else None
        )
        plot_mc_tot(
            ax[0], 
            edges, 
            {'nom': ymc_weight, 'up': ymc_up_weight, 'down': ymc_down_weight}, 
            zorder=-2, 
            label=f'Before {weight}', 
            unc=False, 
            fill=False, 
            color='blue',
            divide=widths if variable_binwidth else None
        )
        plot_mc_tot(
            ax[0], 
            edges, 
            {'nom': ymc, 'up': ymc_up, 'down': ymc_down}, 
            zorder=-1, 
            label=f'After {weight}', 
            unc=False, 
            fill=False, 
            color='red',
            divide=widths if variable_binwidth else None
        )

        # finalize upper panel
        ax[0].set_yscale("log")
        ax[0].legend(
            loc="upper center",
            frameon=True,
            ncols=2,
            framealpha=0.8,
            fontsize=8,
        )
        ax[0].set_ylim(
            max(0.5, np.min(ymc/widths if variable_binwidth else ymc) / 5), 
            np.max(ymc/widths if variable_binwidth else ymc) * 50
        )
        ax[0].set_ylabel("Events/(bin width)" if variable_binwidth else "Events")

        # lower panel
        plot_ratio(
            ax[1], 
            x, 
            edges, 
            [
                {'nom': ymc_weight, 'up': ymc_up_weight, 'down': ymc_down_weight, 'color': 'blue'},
                {'nom': ymc, 'up': ymc_up, 'down': ymc_down, 'color': 'red'}
            ], 
            {'nom': ydata, 'up': ydata_up, 'down': ydata_down}, 
        )

        ax[1].set_ylabel("MC / DATA")
        if variable_label:
            ax[1].set_xlabel(variable_label)
        else:
            ax[1].set_xlabel(variable)

        fig.savefig(
            f"plots/corrections/{region}_{variable}_{weight}.png",
            facecolor="white",
            pad_inches=0.1,
            bbox_inches="tight",
        )
        
        plt.close()


def main():
    analysis_dict = get_analysis_dict()
    samples = analysis_dict["samples"]

    regions = analysis_dict["regions"]
    variables = analysis_dict["variables"]
    nuisances = analysis_dict["nuisances"]
    check_weights = list(analysis_dict["check_weights"].keys())

    colors = analysis_dict["colors"]
    plot_label = analysis_dict["plot_label"]
    year_label = analysis_dict.get("year_label", "Run-II")
    lumi = analysis_dict["lumi"]
    # plot_ylim_ratio = analysis_dict["plot_ylim_ratio"]
    print("Doing plots")

    proc = subprocess.Popen(
        "mkdir -p plots && " + f"cp {get_fw_path()}/data/common/index.php plots/",
        shell=True,
    )
    proc.wait()

    # FIXME add nuisance for stat
    nuisances["stat"] = {
        "name": "stat",
        "type": "stat",
        "samples": dict((skey, "1.00") for skey in samples),
    }

    cpus = 10

    file__ = "histos.root"
    if len(sys.argv) > 1:
        file__ = sys.argv[1]

    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        tasks = []

        input_file = uproot.open(file__)
        for region in regions:
            for variable in variables:
                if "axis" not in variables[variable]:
                    continue
                tasks.append(
                    executor.submit(
                        plot,
                        input_file,
                        region,
                        variable,
                        check_weights,
                        samples,
                        nuisances,
                        lumi,
                        colors,
                        year_label,
                        variables[variable].get("label"),
                    )
                )
        concurrent.futures.wait(tasks)
        for task in tasks:
            task.result()


if __name__ == "__main__":
    main()
