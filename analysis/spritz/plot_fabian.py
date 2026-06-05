import concurrent.futures
import json
import subprocess
import sys
from copy import deepcopy

import matplotlib as mpl
import mplhep as hep
import numpy as np
import math
import uproot
import hist
from spritz.framework.framework import get_analysis_dict, get_fw_path

mpl.use("Agg")
from matplotlib import pyplot as plt

d = deepcopy(hep.style.CMS)

d["font.size"] = 12
d["figure.figsize"] = (5, 5)

d_multidim = deepcopy(d)
d_multidim["font.size"] = 10
d_multidim["figure.figsize"] = (10, 10)

plt.style.use(d)


class Variation(object):

    def __init__(self, name, variations={}, kind=None):
        self.name = name
        self.variations = variations
        self.kind = kind

    def __getitem__(self, key):
        if isinstance(key, slice):
            return Variation(
                name=self.name,
                variations={k: self.variations[k][key.start:key.stop] for k in self.variations},
                kind=self.kind
            )
        elif key in self.variations.keys():
            return self.variations[key]
        else:
            return getattr(self, key)

    def up(self):
        if self.kind == "envelope":
            up = np.max(np.array([variation for variation in self.variations.values()]), axis=0)
            return np.where(up>0., up, np.zeros_like(up))
        elif self.kind == "square":
            varied = np.array([variation for variation in self.variations.values()])
            variation_histo = sum(np.square(varied))
            return np.sqrt(variation_histo)
        elif self.kind == "stdev":
            varied = np.array([variation for variation in self.variations.values()])
            return np.std(varied, axis=0)
        else:
            return self.variations["up"]

    def down(self):
        if self.kind == "envelope":
            down = np.min(np.array([variation for variation in self.variations.values()]), axis=0)
            return np.where(down<0., down, np.zeros_like(down))
        elif self.kind == "square":
            varied = np.array([variation for variation in self.variations.values()])
            variation_histo = sum(np.square(varied))
            return np.sqrt(variation_histo)
        elif self.kind == "stdev":
            varied = np.array([variation for variation in self.variations.values()])
            return np.std(varied, axis=0)
        else:
            return self.variations["down"]

    @classmethod
    def make_variation(cls, directory, nuisance, sample):
        h = directory[f"histo_{sample}"].to_hist()

        name = nuisance["name"]
        type = nuisance["type"]
        kind = nuisance.get("kind")
        samples = nuisance["samples"]
        
        if type == "lnN":
            scaling = float(samples[sample])
            variations = {
                "up": (scaling-1)*h.values(), 
                "down": (1-1.0/scaling)*h.values()
            }

        elif type == "stat":
            variations = {
                "up": np.sqrt(h.variances()), 
                "down": np.sqrt(h.variances())
            }

        elif kind in ["envelope", "square", "stdev"]:
            read_tag = lambda v : v if not isinstance(v,tuple) else v[1]
            variation_tags = [read_tag(v) for v in samples[sample]]
            if sample=="Single Top" and name=="PDFweight":
                variations = {
                    variation: directory[f"histo_{sample}_{variation}"].values().copy() for variation in variation_tags
                }
            else:
                variations = {
                    variation: directory[f"histo_{sample}_{variation}"].values()-h.values() for variation in variation_tags
                }
        
        else:
            variations = {
                "up": abs(directory[f"histo_{sample}_{name}Up"].values()-h.values()),
                "down": abs(directory[f"histo_{sample}_{name}Down"].values()-h.values())
            }
        
        return Variation(name, variations, kind)

    @classmethod
    def add(cls, summands):
        name,kind = summands[0].name,summands[0].kind
        assert all([s.name==name and s.kind==kind for s in summands])

        variations = {v: sum([s.variations[v] for s in summands]) for v in summands[0].variations}
        
        return cls(name=name, variations=variations, kind=kind)



class Histogram(object):

    def __init__(self, name, nom, variations={}, is_data=False, color="black", axis=None):
        self.name = name
        self.nom = nom
        self.variations = variations
        self.is_data = is_data
        self.color = color
        self.axis = axis

    def __getitem__(self, key):
        if isinstance(key, slice):
            edges = self.edges[key.start:key.stop+1]
            if isinstance(self.axis, hist.axis.Variable):
                axis = hist.axis.Variable(edges, name=self.axis.name)
            else:
                axis = hist.axis.Regular(len(edges)-1, edges[0], edges[-1], name=self.axis.name)

            variations = {
                nuis: self.variations[nuis][key.start:key.stop] for nuis in self.variations
            }

            return Histogram(
                name=self.name,
                nom=self.nom[key.start:key.stop],
                variations=variations,
                is_data=self.is_data,
                color=self.color,
                axis=axis
            )
        
        else:
            return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    @property
    def centers(self): return self.axis.centers

    @property
    def x(self): return self.centers

    @property
    def edges(self): return self.axis.edges

    @property
    def widths(self): return self.axis.widths

    @property
    def variable_width(self): return isinstance(self.axis, hist.axis.Variable)

    @property
    def variation_names(self): return list(self.variations.keys())

    def up(self, variations=None):
        if variations is None: variations = self.variation_names
        return np.sqrt(sum([np.square(self.variations[v].up()) for v in variations]))

    def down(self, variations=None):
        if variations is None: variations = self.variation_names
        return np.sqrt(sum([np.square(self.variations[v].down()) for v in variations]))

    def set_axis(self, axis):
        assert len(axis.centers) == len(self.axis.centers)
        self.axis = axis

    def max(self, divide=None, with_unc=False): 
        if divide is None:
            divide = self.widths if self.variable_width else np.ones_like(self.centers)
        if with_unc: 
            return np.max((self.nom+self.up())/divide)
        else: 
            return np.max((self.nom)/divide)
        
    def min(self, divide=None, with_unc=False): 
        if divide is None:
            divide = self.widths if self.variable_width else np.ones_like(self.centers)
        if with_unc: 
            return np.min((self.nom-self.down())/divide)
        else: 
            return np.min((self.nom)/divide)

    @classmethod
    def make_hist(cls, directory, nuisances, sample, is_data=False, color="black"):
        
        h = directory[f"histo_{sample}"].to_hist()
        
        nuisances = {
            k:v for k,v in nuisances.items() if sample in v["samples"] and not v["type"] in ["rateParam","auto"]
        }

        variations = {
            nuis: Variation.make_variation(directory, nuisances[nuis], sample) for nuis in nuisances
        }
        
        return cls(name=sample, nom=h.values().copy(), variations=variations, is_data=is_data, color=color, axis=h.axes[0])

    @classmethod
    def sum_hist(cls, name, histos, is_data=False, color="black"):
        nom = sum([h.nom for h in histos])

        nuisances = union([h.variations.keys() for h in histos])
        variations = {
            nuis: Variation.add(
                [h.variations[nuis] for h in histos if nuis in h.variations]
            ) for nuis in nuisances
        }
        
        return cls(name=name, nom=nom, variations=variations, is_data=is_data, color=color, axis=histos[0].axis)

    @classmethod
    def empty_like(cls, histo, name=None, is_data=None, color=None):
        return cls(
            name=histo.name if name is None else name,
            nom=np.zeros_like(histo.nom),
            variations={},
            is_data=histo.is_data if is_data is None else is_data, 
            color=histo.color if color is None else color, 
            axis=histo.axis
        )

    def plot_data(self, ax, zorder=1, label=False, divide=None):
        if divide is None:
            divide = self.widths if self.variable_width else np.ones_like(self.centers)
        ax.errorbar(
            self.centers,
            self.nom/divide,
            yerr=(self.down()/divide, self.up()/divide),
            fmt="o",
            markersize=4,
            label=f"{self.name} [{int(round(np.sum(self.nom), 0))}]" if label else None, 
            zorder=zorder,
            color=self.color
        )

    def plot_mc_unc(self, ax, zorder=1, label=None, color=None, shaded=False, divide=None, uncertainties=None):
        if divide is None:
            divide = self.widths if self.variable_width else np.ones_like(self.centers)
        unc_up = round(np.sum(self.up(uncertainties)) / np.sum(self.nom) * 100, 2)
        unc_down = round(np.sum(self.down(uncertainties)) / np.sum(self.nom) * 100, 2)
        kwargs = {"alpha": 0.25} if shaded else {"hatch": "///", "facecolor": "none"}
        ax.stairs(
            (self.nom+self.up(uncertainties))/divide,
            edges=self.edges,
            baseline=(self.nom-self.down(uncertainties))/divide,
            label=f"Syst [-{unc_down}, +{unc_up}]%" if label else None,
            fill=True,
            color=self.color if color is None else color, 
            zorder=zorder,
            **kwargs
        )

    def plot_mc(self, ax, baseline=None, zorder=1, draw_unc=False, unc_shaded=False, highlight_unc=[], label=False, label_sum=False, label_unc=False, color=None, fill=False, draw_edge=False, linestyle="solid", divide=None, alpha=1.):
        if divide is None:
            divide = self.widths if self.variable_width else np.ones_like(self.centers)
        
        if draw_unc:
            self.plot_mc_unc(ax, zorder=zorder+(-0.1 if unc_shaded else 0.1), color=darker_color(self.color if color is None else color), shaded=unc_shaded, divide=divide, label=label_unc)

            unc_colors = ["red", "green"]+list(mpl.colors.TABLEAU_COLORS.values())
            for i,unc in enumerate(highlight_unc):
                self.plot_mc_unc(ax, zorder=zorder+(-0.1 if unc_shaded else 0.1), color=unc_colors[i%12], shaded=True, divide=divide, label=label_unc, uncertainties=[unc])
            
        if label:
            if label_sum:
                labeltxt = f"{self.name} [{int(round(np.sum(self.nom), 0))}]"
            else:
                labeltxt = self.name
        else:
            labeltxt = None
        ax.stairs(
            (self.nom if baseline is None else self.nom+baseline)/divide, 
            edges=self.edges, 
            label=labeltxt,
            color=self.color if color is None else color,
            edgecolor=darker_color(self.color if color is None else color),
            linewidth=1. if draw_edge else 0.,
            fill=fill,
            linestyle=linestyle,
            zorder=zorder,
            alpha=alpha
        )

class StackedHistogram(object):

    def __init__(self, name, histos):
        self.name = name
        self.axis = histos[0].axis
        self.histos = []
        
        for h in histos:
            assert h.axis == self.axis
            self.histos.append(h)

    @property
    def centers(self): return self.axis.centers

    @property
    def x(self): return self.centers

    @property
    def edges(self): return self.axis.edges

    @property
    def widths(self): return self.axis.widths

    @property
    def variable_width(self): return isinstance(self.axis, hist.axis.Variable)

    def add(self, histo, position=None):
        assert histo.axis == self.axis
        if position is None:
            self.histos.append(histo)
        else:
            self.histos = self.histos[:position] + [histo] + self.histos[position:]

    def set_axis(self, axis):
        assert len(axis.centers) == len(self.axis.centers)
        for h in self.histos:
            h.set_axis(axis)
        self.axis = axis

    def __getitem__(self, key):
        if isinstance(key, slice):
            return StackedHistogram(
                name=self.name,
                histos=[h[key.start:key.stop] for h in self.histos]
            )
        elif isinstance(key, int):
            return self.histos[key]
        else:
            for h in self.histos:
                if h.name==key:
                    return h
            return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __iter__(self):
        for h in self.histos:
            yield h

    def plot_stack(self, ax, draw_edge=True, zorder=1, label=False, label_sum=True):
        for i,h in enumerate(self.histos):
            h.plot_mc(
                ax, 
                baseline=sum([h["nom"] for h in self.histos[:i]]),
                label=label,
                label_sum=label_sum,
                fill=True,
                color=h.color,
                draw_edge=draw_edge,
                zorder=zorder-i/(len(self.histos)+1)
            )

    def sum(self, name=None, color="black"):
        return Histogram.sum_hist(
            self.name if name is None else name,
            self.histos,
            color=color
        )

    def max(self, divide=None, with_unc=False, total=False):
        if total:
            return max([h.max(divide, with_unc) for h in self.histos])
        else:
            return self.histos[0].max(divide, with_unc)

    def min(self, divide=None, with_unc=False, total=False):
        if total:
            return min([h.min(divide, with_unc) for h in self.histos])
        else:
            return self.histos[0].min(divide, with_unc)


def darker_color(color):
    rgb = list(mpl.colors.to_rgba(color)[:-1])
    darker_factor = 4 / 5
    rgb[0] = rgb[0] * darker_factor
    rgb[1] = rgb[1] * darker_factor
    rgb[2] = rgb[2] * darker_factor
    return tuple(rgb)


def union(lists):
    res = list()
    for l in lists:
        res += [x for x in l if not x in res]
    return res


def main_panel(ax, histos, draw_mc_unc=True, labels=[], labels_unc=[], highlight_unc=[], print_unc=False, label_sum=True, mc_alpha=1.):
    for key,histo in histos.items():
        if isinstance(histo, StackedHistogram):
            histo.plot_stack(ax, label=key in labels, label_sum=label_sum)
        elif histo.is_data:
            histo.plot_data(ax, label=key in labels)
        else:
            histo.plot_mc(ax, draw_unc=draw_mc_unc, label=key in labels, label_unc=key in labels_unc, draw_edge=True, label_sum=label_sum, alpha=mc_alpha)

            if print_unc:
                unc_up = round(np.sum(histo.up()) / np.sum(histo.nom) * 100, 2)
                unc_down = round(np.sum(histo.down()) / np.sum(histo.nom) * 100, 2)
                text = ax.text(
                    0.7, 0.75, fontsize=6, transform=ax.transAxes,
                    s=f"Syst [-{unc_down}, +{unc_up}]%"
                )
                i_highlight = 0
                unc_colors = ["red","green"]+list(mpl.colors.TABLEAU_COLORS.values())
                for syst in histo.variation_names:
                    up = round(np.sum(histo.up([syst])) / np.sum(histo.nom) * 100, 2)
                    down = round(np.sum(histo.down([syst])) / np.sum(histo.nom) * 100, 2)
                    if syst in highlight_unc:
                        text = ax.annotate(
                            text=f"{syst} [-{down}, +{up}]%", xycoords=text, xy=(0,0), 
                            verticalalignment="top", fontsize=6, color=darker_color(unc_colors[i_highlight%12])
                            #verticalalignment="top", fontsize=6, color=darker_color(unc_colors[0])
                        )
                        i_highlight += 1
                    else: 
                        text = ax.annotate(
                            text=f"{syst} [-{down}, +{up}]%", xycoords=text, xy=(0,0), 
                            verticalalignment="top", fontsize=6, color="black"
                        )

def ratio_panel(ax, histos, numerators, denominator, draw_mc_unc=True, printout=False, highlight_unc=[], mc_alpha=1.):
    
    denominator_nom = np.where(histos[denominator]["nom"] >= 1e-6, histos[denominator]["nom"], 1e-6)
    
    # plot denominator
    histos[denominator].plot_mc(ax, draw_unc=draw_mc_unc, unc_shaded=True, highlight_unc=highlight_unc, draw_edge=True, linestyle="dashed", divide=denominator_nom, color="black", alpha=mc_alpha)

    # plot numerators
    if len(numerators)>1:
        offsets = np.linspace(-0.3, 0.3, len(numerators))
    else:
        offsets = np.array([0.0])
    
    res_dict = {}

    for i,numerator in enumerate(numerators):
        edge, x = 0., histos[denominator].centers
        x_i = x.copy()
        for j in range(len(histos[denominator].centers)):
            width = x[j]-edge
            edge = x[j]+width
            x_i[j] = x[j]+offsets[i]*width
        
        if histos[numerator]["is_data"]:
            histos[numerator].plot_data(ax, divide=denominator_nom)
        else:
            histos[numerator].plot_mc(ax, draw_unc=draw_mc_unc, unc_shaded=True, divide=denominator_nom, draw_edge=True, alpha=mc_alpha)
        
        if printout:
            print(json.dumps({
                "edges": list(histos[numerator].edges),
                "nom": list(histos[numerator].nom/denominator_nom),
                "up": list(np.sqrt(np.square(histos[numerator].down())+np.square(histos[denominator].up()))/denominator_nom),
                "down": list(np.sqrt(np.square(histos[numerator].up())+np.square(histos[denominator].down()))/denominator_nom)
            }))


def make_plots(axes, histos, numerators, denominator, draw_mc_unc=True, ncols_legend=3, printout=False, labels=[], labels_unc=[], highlight_unc=[], print_unc=False, xlabel="x", unit=None, legend=True, hide_xlabel=False, hide_ylabel=False, xlog=False, yrange=None, ryrange=None, label_sum=True, mc_alpha=1.):
    main_panel(
        ax=axes[0], 
        histos=histos,
        draw_mc_unc=draw_mc_unc,
        labels=labels, 
        labels_unc=labels_unc,
        highlight_unc=highlight_unc,
        print_unc=print_unc,
        label_sum=label_sum,
        mc_alpha=mc_alpha
    )

    # finalize upper panel
    axes[0].set_yscale("log")
    if legend:
        axes[0].legend(
            loc="upper center",
            frameon=True,
            ncols=ncols_legend,
            framealpha=0.8,
            fontsize=8,
        )

    variable_binwidth = histos[denominator].variable_width
    if yrange is None:
        y0 = max(
            1e-4 if variable_binwidth else 0.5, 
            0.5*min(h.min() for h in histos.values())
        )
        y1 = max(
            h.max() for h in histos.values()
        ) * 1e3
        yrange = (y0, y1)

    ylabel = f"Events / {unit if unit is not None else xlabel}" if variable_binwidth else "Events"
    axes[0].set_ylim(yrange[0], yrange[1])
    axes[0].set_ylabel("" if hide_ylabel else ylabel)
    axes[0].tick_params(labelbottom=False)
    axes[0].set_xlabel("")

    # lower panel
    ratio_panel(
        axes[1], 
        histos=histos,
        numerators=numerators,
        denominator=denominator,
        draw_mc_unc=draw_mc_unc,
        printout=printout,
        highlight_unc=highlight_unc,
        mc_alpha=mc_alpha
    )

    if ryrange is None:
        denominator_nom = np.where(histos[denominator]["nom"] >= 1e-6, histos[denominator]["nom"], 1e-6)
        ry0 = max(0., 
            min([h.min(divide=denominator_nom) for h in histos.values()]
                + [histos[denominator].min(divide=denominator_nom, with_unc=True)]))
        ry1 = min(2.,
            max([h.max(divide=denominator_nom) for h in histos.values()]
                + [histos[denominator].max(divide=denominator_nom, with_unc=True)]))
        rydiff = 1.1*max(ry1-1., 1.-ry0)
        ryrange = (1-rydiff, 1+rydiff)

    xlabel = xlabel + f" ({unit})" if unit is not None else xlabel
    axes[1].set_xlim(histos[list(histos.keys())[0]].edges[0], histos[list(histos.keys())[0]].edges[-1])
    axes[1].set_ylim(ryrange[0], ryrange[1])
    axes[1].set_ylabel("" if hide_ylabel else f"Ratio to {denominator}")
    axes[1].set_xlabel("" if hide_xlabel else xlabel)

    if xlog: 
        axes[1].set_xscale("log")
        if axes[1].get_xlim()[0]==0:
            xmin = histos[list(histos.keys())[0]].edges[1]/2
            axes[1].set_xlim(xmin, histos[list(histos.keys())[0]].edges[-1])


def make_plots_multidim(axes, histos, numerators, denominator, h_axis, variable_label, unit, draw_mc_unc=True, highlight_unc=[], xlog=False, ryrange=(0.85,1.15), label_sum=True):
    if variable_label is None:
        variable_label = [h_axis[i].name for i in range(len(h_axis))]

    nbins = len(h_axis[0].centers)
    if len(h_axis)==3:
        ncols = len(h_axis[1].centers)
        nrows = len(h_axis[2].centers)
    elif len(h_axis)==2:
        nrows = math.floor(math.sqrt(len(h_axis[1].centers)))
        ncols = math.ceil(len(h_axis[1].centers)/nrows)

    npanels = int(len(histos[denominator].centers)/nbins)
    widths = np.tile(h_axis[0].widths, npanels)
    y0 = max(
        1e-8, 0.5*min(h.min(divide=widths) for h in histos.values())
    )
    y1 = max(
        h.max(divide=widths) for h in histos.values()
    ) * 100

    bbox = {"boxstyle":"square", "alpha":1.0, "fc":"white", "ec":"black"} 
    finished = False
    for irow in range(nrows):
        for icol in range(ncols):
            histos_sliced = {}
            for key,histo in histos.items():
                if (ncols*nbins)*irow+nbins*(icol+1) > len(histo.widths):
                    finished = True
                    break
                histos_sliced[key] = histo[(ncols*nbins)*irow+nbins*icol:(ncols*nbins)*irow+nbins*(icol+1)]
                histos_sliced[key].set_axis(h_axis[0])
            if finished:
                break

            make_plots(
                (axes[2*irow,icol],axes[2*irow+1,icol]),
                histos=histos_sliced,
                numerators=numerators,
                denominator=denominator,
                draw_mc_unc=draw_mc_unc,
                xlabel=variable_label[0],
                unit=unit[0] if unit is not None else None,
                legend=False,
                hide_xlabel=irow!=nrows-1,
                hide_ylabel=icol!=0,
                highlight_unc=highlight_unc,
                xlog=xlog,
                yrange=(y0,y1),
                ryrange=ryrange,
                label_sum=label_sum
            )

            if len(h_axis)==3:
                axes[2*irow,icol].text(0.96, 0.95,
                    f"${h_axis[1].edges[icol]} <${variable_label[1]}<$ {h_axis[1].edges[icol+1]}$",
                    fontsize=7,
                    horizontalalignment="right",
                    verticalalignment="top",
                    transform=axes[2*irow,icol].transAxes,
                    bbox=bbox
                )

                axes[2*irow,icol].text(0.04, 0.95,
                    f"${h_axis[2].edges[irow]} <${variable_label[2]}$< {h_axis[2].edges[irow+1]}$",
                    fontsize=7,
                    horizontalalignment="left",
                    verticalalignment="top",
                    transform=axes[2*irow,icol].transAxes,
                    bbox=bbox
                )

            elif len(h_axis)==2:
                axes[2*irow,icol].text(0.96, 0.95,
                    f"${h_axis[1].edges[nrows*irow+icol]} <${variable_label[1]}<$ {h_axis[1].edges[nrows*irow+icol+1]}$",
                    fontsize=7,
                    horizontalalignment="right",
                    verticalalignment="top",
                    transform=axes[2*irow,icol].transAxes,
                    bbox=bbox
                )


def setup_multifig(ncols, nrows):
    fig = plt.figure(dpi=200) 
    ax = np.empty((2*nrows,ncols), dtype=plt.Axes)
    gs = mpl.gridspec.GridSpec(
        3*nrows-1, 2*ncols-1,
        figure=fig,
        height_ratios=([2,1,0.05]*nrows)[:-1],
        width_ratios=([1,0]*ncols)[:-1]
    )
    
    for i in range(nrows):
        for j in range(ncols): 
            ax[2*i,j] = fig.add_subplot(gs[3*i,2*j], sharex=ax[0,0], sharey=ax[0,0])
            ax[2*i+1,j] = fig.add_subplot(gs[3*i+1,2*j], sharex=ax[0,0], sharey=ax[1,0])

    for axij in ax.flat:
        axij.label_outer()

    return fig,ax


def plot(
    input_file,
    region,
    variable,
    samples,
    nuisances,
    lumi,
    colors,
    year_label,
    variable_label=None,
    unit=None,
    xlog=False,
    axis=None,
    plotFakes=False,
    plotVariations=False
):
    print("Doing ", region, variable)

    directory = input_file[f"{region}/{variable}"]
    mc_samples = [x for x in samples if not samples[x].get("is_data", False) and not x in ["W+Jets","QCD","TTToSemiLeptonic"]]
    mcfakes_samples = [x for x in samples if x in ["W+Jets","QCD","TTToSemiLeptonic"]]
    
    # get the histograms
    histos = {
        sample: Histogram.make_hist(directory, nuisances, sample, is_data=samples[sample].get("is_data", False), color=colors.get(sample,"black"))
        for sample in samples
    }

    # prepare total MC histogram
    stack_mc = StackedHistogram("", [histos[sample] for sample in mc_samples])
    histo_mc = stack_mc.sum("Tot MC", color="grey")

    if plotFakes and not "_ss" in region:
        directory_ss = input_file[f"{region}_ss/{variable}"]
        histos_ss = {
            sample: Histogram.make_hist(directory_ss, nuisances, sample, is_data=samples[sample].get("is_data", False), color=colors.get(sample,"black"))
            for sample in samples
        }
        stack_mc_ss = StackedHistogram("", [histos_ss[sample] for sample in mc_samples])
        histo_mc_ss = stack_mc_ss.sum("Tot MC SS", color="grey")

        if "Data" in histos_ss:
            histo_data_ss = histos_ss["Data"]
        else:
            histo_data_ss = Histogram.empty_like(histo_mc_ss, name="Data", is_data=True, color="black")

        # subtract MC from data to get fakes
        variations_fakes = {}
        variations_fakes["stat"] = Variation(
            "stat", {
                "up": np.sqrt(np.square(histo_data_ss.up(["stat"]))+np.square(histo_mc_ss.down(["stat"]))),
                "down": np.sqrt(np.square(histo_data_ss.down(["stat"]))+np.square(histo_mc_ss.up(["stat"]))) }
        )

        if "fakerw_param" in histo_data_ss.variation_names:
            variations_fakes["fakerw_param"] = Variation(
                "fakerw_param", {   
                    "up": histo_data_ss.variations["fakerw_param"].variations["up"]-histo_mc_ss.variations["fakerw_param"].variations["up"],
                    "down": histo_data_ss.variations["fakerw_param"].variations["down"]-histo_mc_ss.variations["fakerw_param"].variations["down"] }, 
                "weight"
            )
        
        if "fakerw_model" in histo_data_ss.variation_names:
            variations_fakes["fakerw_model"] = Variation(
                "fakerw_model", {
                    "fakerw_model_exp": histo_data_ss.variations["fakerw_model"].variations["fakerw_model_exp"]-histo_mc_ss.variations["fakerw_model"].variations["fakerw_model_exp"],
                    "fakerw_model_erf": histo_data_ss.variations["fakerw_model"].variations["fakerw_model_erf"]-histo_mc_ss.variations["fakerw_model"].variations["fakerw_model_erf"] },
                "envelope"
            )

        histo_fakes = Histogram(
            "Fakes", 
            histo_data_ss.nom-histo_mc_ss.nom, 
            variations=variations_fakes, 
            is_data=True, 
            color=colors["Fakes"],
            axis=histo_data_ss.axis
        )

        stack_mc.add(histo_fakes, position=0)
        histo_mc = stack_mc.sum("Tot MC", color="grey")

    # prepare data histogram
    if "Data" in histos:
        histo_data = histos["Data"]
    else:
        histo_data = Histogram.empty_like(histo_mc, name="Data", is_data=True, color="black")

    # make plots
    if isinstance(axis, list):
        plt.style.use(d_multidim)
        if len(axis)==3:
            ncols = len(axis[1].centers)
            nrows = len(axis[2].centers)
        elif len(axis)==2:
            nrows = math.floor(math.sqrt(len(axis[1].centers)))
            ncols = math.ceil(len(axis[1].centers)/nrows)

        fig, ax = setup_multifig(ncols, nrows)
        fig.tight_layout(pad=-0.4)
        hep.cms.label(region, rlabel="", data=True, ax=ax[0,0])
        hep.label.exp_label(data=True, lumi=round(lumi, 2), year=year_label, ax=ax[0,-1])

        make_plots_multidim(
            ax, 
            histos={"MC Stack": stack_mc, "MC": histo_mc, "Data": histo_data},
            numerators=["Data"],
            denominator="MC", 
            h_axis=axis, 
            variable_label=variable_label, 
            unit=unit,
            highlight_unc=[],
            xlog=xlog,
            ryrange=(0.8,1.2)
        )
    else:
        plt.style.use(d)
        fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3,1]}, dpi=200)
        hep.cms.label(region, data=True, lumi=round(lumi, 2), ax=ax[0], year=year_label)
        fig.tight_layout(pad=-0.5)
        make_plots(
            ax,
            histos={"MC Stack": stack_mc, "MC": histo_mc, "Data": histo_data},
            numerators=["Data"],
            denominator="MC",
            labels=["MC Stack","MC","Data"],
            labels_unc=["MC"],
            xlabel=variable_label if variable_label is not None else variable,
            unit=unit,
            highlight_unc=[],
            print_unc=False,
            xlog=xlog,
            ryrange=(0.8,1.2)
        )
    
    fig.savefig(
        f"plots/{region}_{variable}.pdf",
        facecolor="white",
        pad_inches=0.1,
        bbox_inches="tight",
    )

    plt.close()
    

    if plotVariations and variable in ["mll_medium"]:#,"triple_diff"]:

        for variation in histo_mc.variation_names:
            name = nuisances[variation]["name"]
            type = nuisances[variation]["type"]
            kind = nuisances[variation].get("kind")
            samples = nuisances[variation]["samples"]

            histo_mc_nominal = Histogram.empty_like(histo_mc, name=f"{name}_nominal", color="black")
            histo_mc_nominal.nom = histo_mc.nom
            histo_mc_nominal.variations = { variation: histo_mc.variations[variation] }
            
            histo_mc_variation = {}
            var_colors = ["red","blue","green","purple","cyan","magenta","grey","brown","orange"]
            if kind in ["envelope","square","stdev"]:
                for i,key in enumerate(histo_mc.variations[variation].variations.keys()):
                    histo_mc_variation[key] = Histogram.empty_like(histo_mc, name=key)
                    histo_mc_variation[key].nom = histo_mc.nom + histo_mc.variations[variation].variations[key]
                    histo_mc_variation[key].color = var_colors[i % 9]
            else:
                histo_mc_variation["up"] = Histogram.empty_like(histo_mc, name=f"{name}_up")
                histo_mc_variation["down"] = Histogram.empty_like(histo_mc, name=f"{name}_down")
                histo_mc_variation["up"].nom = histo_mc.nom + histo_mc.variations[variation].up()
                histo_mc_variation["down"].nom = histo_mc.nom - histo_mc.variations[variation].down()
                histo_mc_variation["up"].color = var_colors[0]
                histo_mc_variation["down"].color = var_colors[1]

            if isinstance(axis, list):
                plt.style.use(d_multidim)
                if len(axis)==3:
                    ncols = len(axis[1].centers)
                    nrows = len(axis[2].centers)
                elif len(axis)==2:
                    nrows = math.floor(math.sqrt(len(axis[1].centers)))
                    ncols = math.ceil(len(axis[1].centers)/nrows)

                fig, ax = setup_multifig(ncols, nrows)
                fig.tight_layout(pad=-0.4)
                hep.cms.label(region, rlabel="", data=True, ax=ax[0,0])
                hep.label.exp_label(data=True, lumi=round(lumi, 2), year=year_label, ax=ax[0,-1])

                make_plots_multidim(
                    ax, 
                    histos={f"{name}_nominal": histo_mc_nominal} | {k: histo_mc_variation[k] for k in histo_mc.variations[variation].variations.keys()},
                    numerators=[k for k in histo_mc.variations[variation].variations.keys()],
                    denominator=f"{name}_nominal", 
                    h_axis=axis, 
                    variable_label=variable_label, 
                    unit=unit,
                    highlight_unc=[],
                    xlog=xlog,
                    ryrange=(0.8,1.2),
                    label_sum=False
                )
            else:
                plt.style.use(d)
                fig, ax = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2,3]}, dpi=200)
                hep.cms.label(region, data=True, lumi=round(lumi, 2), ax=ax[0], year=year_label)
                fig.tight_layout(pad=-0.5)
                
                make_plots(
                    ax,
                    histos={f"{name}_nominal": histo_mc_nominal} | {k: histo_mc_variation[k] for k in histo_mc_variation.keys()},
                    numerators=[k for k in histo_mc_variation.keys()],
                    denominator=f"{name}_nominal",
                    labels=[f"{name}_nominal"]+[k for i,k in enumerate(histo_mc_variation.keys()) if i<12],
                    xlabel=variable_label if variable_label is not None else variable,
                    unit=unit,
                    xlog=xlog,
                    draw_mc_unc=True,
                    label_sum=False,
                    ncols_legend=4,
                    mc_alpha=0.3 if len(histo_mc_variation.keys())>10 else 1.
                )

            fig.savefig(
                f"plots/variations/{region}_{variable}_{name}.pdf",
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

    colors = analysis_dict["colors"]
    plot_label = analysis_dict["plot_label"]
    year_label = analysis_dict.get("year_label", "Run-II")
    lumi = analysis_dict["lumi"]

    # FIXME add nuisance for stat
    nuisances["stat"] = {
        "name": "stat",
        "type": "stat",
        "samples": dict((skey, "1.00") for skey in samples),
    }

    plotFakes = (len(sys.argv)>1 and sys.argv[1]=="--fakes") or (len(sys.argv)>2 and sys.argv[2]=="--fakes")
    plotVariations = (len(sys.argv)>1 and sys.argv[1]=="--variations") or (len(sys.argv)>2 and sys.argv[2]=="--variations")

    cmd_mkdir = f"mkdir -p plots && cp {get_fw_path()}/data/common/index.php plots/"
    if plotVariations:
        cmd_mkdir += f" && mkdir -p plots/variations && cp {get_fw_path()}/data/common/index.php plots/variations/"
    
    proc = subprocess.Popen(cmd_mkdir, shell=True)
    proc.wait()

    cpus = 15

    print("Doing plots")

    with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
        tasks = []

        input_file = uproot.open("histos.root")
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
                        samples,
                        nuisances,
                        lumi,
                        colors,
                        year_label,
                        variables[variable].get("label"),
                        variables[variable].get("unit"),
                        variables[variable].get("xlog", False),
                        variables[variable].get("axis"),
                        plotFakes,
                        plotVariations
                    )
                )
        concurrent.futures.wait(tasks)
        for task in tasks:
            task.result()


if __name__ == "__main__":
    main()
