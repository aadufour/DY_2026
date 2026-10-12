import gc
import io
import json
import os
import sys
import time
import traceback as tb
import zlib
from copy import deepcopy
from itertools import combinations

import awkward as ak
import cloudpickle
import hist
import numpy as np
import uproot

from matplotlib.colors import LinearSegmentedColormap, to_hex
from spritz.lookup_tools.eft_operator_indices import get_group, eft_rw_indices


def get_fw_path():
    path_fw = os.getenv("SPRITZ_PATH")
    if path_fw is None:
        raise Exception("Could not find SPRITZ_PATH variable, remember to source!")
    return path_fw


get_fw_path()


def get_config_path():
    path = os.path.abspath(".")
    print("Working in analysis path:", path)
    return path


def get_analysis_dict(path=None):
    if not path:
        path = get_config_path()
    sys.path.insert(0, path)

    exec("import config as analysis_cfg", globals(), globals())

    return analysis_cfg.__dict__  # type: ignore # noqa: F821


def get_batch_cfg():
    # batch_config.local.json (not tracked) overrides batch_config.json key
    # by key, for site-specific settings
    batch_cfg = dict()
    for name in ["batch_config.json", "batch_config.local.json"]:
        if os.path.isfile(f"{get_fw_path()}/{name}"):
            with open(f"{get_fw_path()}/{name}", "r") as file:
                batch_cfg.update(json.load(file))
    return {
        "X509_USER_PROXY": batch_cfg.get("X509_USER_PROXY", None),
        # false: don't let condor transfer the proxy (the job only exports its
        # path), for pools where proxy delegation fails and the path is shared
        "X509_DELEGATE": batch_cfg.get("X509_DELEGATE", True),
        "SINGULARITY_IMAGE": batch_cfg.get("SINGULARITY_IMAGE", None),
        "BATCH_SYSTEM": batch_cfg.get("BATCH_SYSTEM", "condor"),
        "JOB_FLAVOUR": batch_cfg.get("JOB_FLAVOUR", None),
        "MACHINES": batch_cfg.get("MACHINES", []),
        "REQUEST_MEMORY": batch_cfg.get("REQUEST_MEMORY", 2048),
        # memory for jobs with EFT chunks (None: same as REQUEST_MEMORY)
        "REQUEST_MEMORY_EFT": batch_cfg.get("REQUEST_MEMORY_EFT", None),
        # stream out.txt/err.txt to the submit node while the job runs (condor)
        "STREAM_OUTPUT": batch_cfg.get("STREAM_OUTPUT", False),
    }


def correctionlib_wrapper(ceval):
    return ceval.evaluate


def max_vec(vec, val):
    return ak.where(vec > val, vec, val)


def over_under(val, min, max):
    val = ak.where(val >= max, max, val)
    val = ak.where(val <= min, min, val)
    return val


def m_pi_pi(phi):
    return ak.where(
        phi > np.pi,
        phi - 2 * np.pi,
        ak.where(
            phi <= -np.pi,
            phi + 2 * np.pi,
            phi,
        ),
    )


def read_events(filename, start=0, stop=100, read_form={}):
    print("start reading", flush=True)
    uproot_options = dict(
        timeout=30,
        handler=uproot.source.xrootd.XRootDSource,
        num_workers=1,
        use_threads=False,
    )
    f = uproot.open(filename, **uproot_options)
    tree = f["Events"]
    start = min(start, tree.num_entries)
    stop = min(stop, tree.num_entries)
    if start >= stop:
        return ak.Array([])

    branches = [k.name for k in tree.branches]

    events = {}
    form = deepcopy(read_form)

    all_branches = []
    for coll in form:
        coll_branches = form[coll]["branches"]
        if len(coll_branches) == 0:
            if coll in branches:
                all_branches.append(coll)
        else:
            for branch in coll_branches:
                branch_name = coll + "_" + branch
                if branch_name in branches:
                    all_branches.append(branch_name)

    events_bad_form = tree.arrays(
        all_branches,
        entry_start=start,
        entry_stop=stop,
        decompression_executor=uproot.source.futures.TrivialExecutor(),
        interpretation_executor=uproot.source.futures.TrivialExecutor(),
    )
    f.close()

    for coll in form:
        d = {}
        coll_branches = form[coll].pop("branches")

        if len(coll_branches) == 0:
            if coll in branches:
                events[coll] = events_bad_form[coll]
            continue

        for branch in coll_branches:
            branch_name = coll + "_" + branch
            if branch_name in branches:
                if branch_name.endswith("phi"):
                    vals = events_bad_form[branch_name]
                    vals = over_under(vals, -np.pi, np.pi)
                    d[branch] = vals
                else:
                    d[branch] = events_bad_form[branch_name]

        if len(d.keys()) == 0:
            print("did not find anything for", coll, filename, file=sys.stderr)
            continue

        events[coll] = ak.zip(d, **form[coll])
        del d

    _events = ak.zip(events, depth_limit=1)
    del events
    gc.collect()
    return _events


def add_dict(d1, d2):
    if isinstance(d1, dict):
        d = {}
        common_keys = d1.keys() & d2.keys()
        for key in common_keys:
            if key in ("eft_names", "eft_batch_size"):
                # identical across every chunk of the same dataset
                d[key] = d1[key]
            else:
                d[key] = add_dict(d1[key], d2[key])
        for key in d1.keys()-common_keys:
            d[key] = d1[key]
        for key in d2.keys()-common_keys:
            d[key] = d2[key]
        return d
    elif isinstance(d1, ak.Array):
        return ak.concatenate([d1, d2])
    elif isinstance(d1, np.ndarray) and d1.ndim != 0:
        print("Debug np", d1, d2)
        return np.concatenate([d1, d2])
    elif isinstance(d1, set):
        return d1 | d2
    elif isinstance(d1, list):
        # A list of per-batch hist.Hist objects (the megahisto eft_reweighting layout)
        return [add_dict(a, b) for a, b in zip(d1, d2)]
    elif isinstance(d1, EFTHistogram):
        # in place: these can be ~1GB, and every caller (runner in-job merge,
        # spritz-merge reduction) discards d1 afterwards anyway
        d1 += d2
        return d1
    else:
        try:
            return d1 + d2
        except:
            print()
            print('d1')
            print(d1)
            print()
            print('d2')
            print(d2)
            print()


def add_dict_iterable(iterable):
    tmp = None
    for it in iterable:
        if tmp is None:
            tmp = it
        else:
            tmp = add_dict(tmp, it)
    return tmp


class EFTHistogram:
    """Compact accumulator for many EFT histograms sharing one binning: the
    N reweighting points, or the N covariance terms between them.

    Only the contents of real categories are stored, as plain arrays:
        values:    (n_names, n_systs, n_regions, n_bins)
        variances: (n_names, n_regions, n_bins), nominal only (or None)
    where n_bins runs over the flattened variable axes *including* their
    under/overflow bins. to_hist() rebuilds the exact per-name histogram the
    runner used to book: (*axes, category, syst) with Weight storage.
    """

    def __init__(self, axes, regions, systs, n_names, variances=False):
        self.axes = list(axes)
        self.regions = list(regions)
        self.systs = list(systs)
        self.shape = tuple(ax.extent for ax in self.axes)
        n_bins = int(np.prod(self.shape))
        self.values = np.zeros((n_names, len(self.systs), len(self.regions), n_bins))
        self.variances = np.zeros((n_names, len(self.regions), n_bins)) if variances else None

    def bin_index(self, vals):
        """Flat (flow-inclusive) bin index of each event, and a mask of
        events landing in an existing bin (same as hist.fill drops)."""
        flow_idx = []
        valid = np.ones(len(vals[self.axes[0].name]), dtype=bool)
        for ax in self.axes:
            i = np.asarray(ax.index(vals[ax.name]), dtype=np.intp)
            if ax.traits.underflow:
                i = i + 1
            valid &= (i >= 0) & (i < ax.extent)
            flow_idx.append(i)
        flow_idx = [np.clip(i, 0, ax.extent - 1) for i, ax in zip(flow_idx, self.axes)]
        return np.ravel_multi_index(flow_idx, self.shape), valid

    @staticmethod
    def _segments(flat, valid):
        """Order events by bin: returns (order, bins, starts) such that
        order[starts[k]:starts[k+1]] are the events in bin bins[k]."""
        idx = np.flatnonzero(valid)
        order = idx[np.argsort(flat[idx], kind="stable")]
        bins, starts = np.unique(flat[order], return_index=True)
        return order, bins, starts

    def fill(self, vals, weights, region, syst):
        """weights: (n_events, n_names), one column per name."""
        flat, valid = self.bin_index(vals)
        if not valid.any():
            return
        order, bins, starts = self._segments(flat, valid)
        w = weights[order]
        r, s = self.regions.index(region), self.systs.index(syst)
        self.values[:, s, r, bins] += np.add.reduceat(w, starts, axis=0).T
        if self.variances is not None and syst == "nom":
            self.variances[:, r, bins] += np.add.reduceat(w**2, starts, axis=0).T

    def fill_covariance(self, vals, x, idx_i, idx_j, region, syst="nom"):
        """Covariance terms sum_e x_ei * x_ej for the pairs (idx_i, idx_j)
        of columns of x (n_events, n_points), via one Gram matrix per bin
        instead of materializing an (n_events, n_pairs) weight matrix."""
        flat, valid = self.bin_index(vals)
        if not valid.any():
            return
        order, bins, starts = self._segments(flat, valid)
        xs = x[order]
        r, s = self.regions.index(region), self.systs.index(syst)
        stops = np.append(starts[1:], len(order))
        for b, start, stop in zip(bins, starts, stops):
            xb = xs[start:stop]
            gram = xb.T @ xb
            self.values[:, s, r, b] += gram[idx_i, idx_j]

    def __iadd__(self, other):
        if not (
            self.values.shape == other.values.shape
            and self.regions == other.regions
            and self.systs == other.systs
            and self.axes == other.axes
        ):
            raise ValueError("Cannot add EFTHistograms with different layouts")
        self.values += other.values
        if self.variances is not None:
            self.variances += other.variances
        return self

    def to_hist(self, idx):
        h = hist.Hist(
            *self.axes,
            hist.axis.StrCategory(self.regions, name="category"),
            hist.axis.StrCategory(self.systs, name="syst"),
            hist.storage.Weight(),
        )
        n_reg, n_syst = len(self.regions), len(self.systs)
        view = h.view(flow=True)
        # (n_systs, n_regions, n_bins) -> (*axes, n_regions, n_systs)
        values = self.values[idx].reshape(n_syst, n_reg, *self.shape)
        view.value[..., :n_reg, :n_syst] = np.moveaxis(values, (0, 1), (-1, -2))
        if self.variances is not None:
            variances = self.variances[idx].reshape(n_reg, *self.shape)
            view.variance[..., :n_reg, self.systs.index("nom")] = np.moveaxis(variances, 0, -1)
        return h


def expand_eft_combined(results, part=None):
    """Expand EFT "megahisto", where a dataset's N_weights or
    N_covariance histograms are stored as one entry.

    part=(i, n) expands only every n-th EFT histogram starting from the i-th
    (and the non-EFT entries only for i == 0), so that the n parts together
    are the full expansion"""
    expanded = {}
    for key, entry in results.items():
        if key.endswith("__eft_combined"):
            dataset = key[: -len("__eft_combined")]
        elif key.endswith("__eft_points"):
            dataset = key[: -len("__eft_points")]
        elif key.endswith("__eft_covariances"):
            dataset = key[: -len("__eft_covariances")]
        else:
            if part is None or part[0] == 0:
                expanded[key] = entry
            continue
        eft_names = entry["eft_names"]
        batch_size = entry.get("eft_batch_size", 1)  # unused for EFTHistogram
        for idx, name in enumerate(eft_names):
            if part is not None and idx % part[1] != part[0]:
                continue
            batch_idx = idx // batch_size
            local_idx = idx % batch_size
            histos = {}
            for variable, batch_histos in entry["histos"].items():
                if isinstance(batch_histos, EFTHistogram):
                    histos[variable] = batch_histos.to_hist(idx)
                    continue
                batch_histo = batch_histos[batch_idx]
                axis_names = [ax.name for ax in batch_histo.axes]
                sub_pos = axis_names.index("subsample")
                slicer = [slice(None)] * len(batch_histo.axes)
                slicer[sub_pos] = hist.loc(local_idx)
                histos[variable] = batch_histo[tuple(slicer)]
            expanded[f"{dataset}_{name}"] = {
                "sumw": entry["sumw"],
                "nevents": entry["nevents"],
                "events": entry.get("events", 0),
                "histos": histos,
            }
    return expanded


def big_process(process, filenames, start, stop, read_form, **kwargs):
    t0 = time.time()
    events, error = 0, ""
    print(filenames)
    
    for filename in filenames:
        try:
            events = read_events(filename, start=start, stop=stop, read_form=read_form)
            break
        except Exception as e:
            error += "".join(tb.format_exception(None, e, e.__traceback__))
            continue

    if isinstance(events, int):
        print(error, file=sys.stderr)
        raise Exception("Error, could not read any of the filenames\n" + error, filenames)

    t_reading = time.time() - t0
    print(f"created events ({t_reading:.2f} s)", flush=True)
    
    if len(events) == 0:
        return {}
    
    results = {"real_results": 0, "performance": {}}
    results["real_results"] = process(events, **kwargs)
    t_reading = time.time() - t0
    print(f"processed events ({t_reading:.2f} s)", flush=True)
    results["performance"][f"{filename}_{start}"] = {"total": time.time()-t0, "read": t_reading}
    
    del events
    gc.collect()
    return results


class _ZlibReader(io.RawIOBase):
    """Decompress a zlib stream on the fly, so unpickling never holds the
    compressed and decompressed bytes of the whole file in memory"""

    def __init__(self, file, bufsize=1 << 20):
        self._file = file
        self._bufsize = bufsize
        self._d = zlib.decompressobj()

    def readable(self):
        return True

    def readinto(self, b):
        n = len(b)
        if n == 0:
            return 0
        while True:
            if self._d.unconsumed_tail:
                out = self._d.decompress(self._d.unconsumed_tail, n)
            elif self._d.eof:
                return 0
            else:
                data = self._file.read(self._bufsize)
                if not data:
                    raise EOFError("truncated zlib stream")
                out = self._d.decompress(data, n)
            if out:
                b[: len(out)] = out
                return len(out)


class _ZlibWriter(io.RawIOBase):
    """Compress on the fly while pickling, see _ZlibReader"""

    def __init__(self, file):
        self._file = file
        self._c = zlib.compressobj()

    def writable(self):
        return True

    def write(self, b):
        self._file.write(self._c.compress(b))
        return memoryview(b).nbytes

    def finish(self):
        self._file.write(self._c.flush())


def low_memory_executor(max_workers, preload=()):
    """ProcessPoolExecutor whose workers exit after every task, so the memory
    a task freed goes back to the OS instead of staying in the worker's heap.
    That needs a non-fork start method: forkserver, with `preload` modules
    imported once in the server so each new worker starts in ~20 ms"""
    import concurrent.futures
    import multiprocessing

    if sys.version_info < (3, 11):
        return concurrent.futures.ProcessPoolExecutor(max_workers=max_workers)
    ctx = multiprocessing.get_context("forkserver")
    ctx.set_forkserver_preload(list(preload))
    return concurrent.futures.ProcessPoolExecutor(
        max_workers=max_workers, mp_context=ctx, max_tasks_per_child=1
    )


def read_chunks(filename, readable=False, low_memory=False):
    if not readable:
        with open(filename, "rb") as file:
            if low_memory:
                chunks = cloudpickle.load(io.BufferedReader(_ZlibReader(file)))
            else:
                chunks = cloudpickle.loads(zlib.decompress(file.read()))
        return chunks
    else:
        with open(filename, "r") as file:
            chunks = json.load(file)
        return chunks


def write_chunks(d, filename, readable=False, low_memory=False):
    # low_memory streams (de)compression instead of materializing the whole
    # pickle; the file format is the same either way
    if not readable:
        with open(filename, "wb") as file:
            if low_memory:
                writer = _ZlibWriter(file)
                cloudpickle.dump(d, writer)
                writer.finish()
            else:
                file.write(zlib.compress(cloudpickle.dumps(d)))
    else:
        with open(filename, "w") as file:
            json.dump(d, file)


def get_rw_idx(dataset, point):
    group = get_group(dataset)
    rw_indices = eft_rw_indices[group]
    return rw_indices[point]


def get_eft_points(operators, linear=False):
    eft_points = ["sm"]
    eft_points += [f"w1_{op}" for op in operators]
    eft_points += [f"wm1_{op}" for op in operators]
    if not linear:
        eft_points += [f"w11_{op1}_{op2}" for op1, op2 in combinations(operators, 2)]
    return eft_points


def interpolate_colors(base_colors, n_colors):
    """
    Interpolate a list of hex colors.

    Parameters
    ----------
    base_colors : list[str]
        List of hex colors, e.g. ["#ff0000", "#00ff00"]
    n_colors : int
        Number of output colors requested

    Returns
    -------
    list[str]
        Interpolated hex colors
    """

    cmap = LinearSegmentedColormap.from_list(
        "custom_cmap",
        base_colors,
        N=n_colors,
    )

    return [to_hex(cmap(i / (n_colors - 1))) for i in range(n_colors)]


# plots
cmap_petroff = [ # https://github.com/mpetroff/accessible-color-cycles
    "#3f90da",
    "#ffa90e",
    "#bd1f01",
    "#94a4a2",
    "#832db6",
    "#a96b59",
    "#e76300",
    "#b9ac70",
    "#717581",
    "#92dadd"
]
cmap_pastel = [
    "#A1C9F4",
    "#FFB482",
    "#8DE5A1",
    "#FF9F9B",
    "#D0BBFF",
    "#DEBB9B",
    "#FAB0E4",
    "#CFCFCF",
    "#FFFEA3",
    "#B9F2F0",
]
