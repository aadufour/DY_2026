import pickle, sys
import numpy as np

OP = sys.argv[1] if len(sys.argv) > 1 else 'cHDD'

with open('/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr.pkl', 'rb') as f:
    prop = pickle.load(f)
with open('/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_parallel.pkl', 'rb') as f:
    base = pickle.load(f)

edges = np.arange(80, 102, 0.5)

def shape_and_n(cache, op):
    w = cache['w_p1'][op]
    mll = cache['mll']
    h, _ = np.histogram(mll, bins=edges, weights=w)
    n, _ = np.histogram(mll, bins=edges)
    return h / h.sum(), n

sp, n_prop = shape_and_n(prop, OP)
sb, n_base = shape_and_n(base, OP)

print(f"Operator: {OP} (Wilson coefficient = +1, all others = 0)")
print(f"{'mll bin':>14} {'propcorr':>10} {'baseline':>10} {'ratio':>8} {'N_prop':>8} {'N_base':>8}")
for i in range(len(edges) - 1):
    lo, hi = edges[i], edges[i + 1]
    r = sp[i] / sb[i] if sb[i] > 0 else float('nan')
    print(f"{lo:6.1f}-{hi:<6.1f} {sp[i]:10.5f} {sb[i]:10.5f} {r:8.4f} {n_prop[i]:8d} {n_base[i]:8d}")
