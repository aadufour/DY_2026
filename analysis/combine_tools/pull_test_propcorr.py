import pickle, sys
import numpy as np

OP = sys.argv[1] if len(sys.argv) > 1 else 'cHDD'

with open('/grid_mnt/data__data.polcms/cms/adufour/LHE/propcorr/CACHE/lhe_cache_propcorr.pkl', 'rb') as f:
    prop = pickle.load(f)
with open('/grid_mnt/data__data.polcms/cms/adufour/LHE/SYST_slc7/CACHE/lhe_cache_parallel.pkl', 'rb') as f:
    base = pickle.load(f)

edges = np.arange(80, 102, 0.5)

def shape_and_n(cache, op):
    w = cache['w_p1'][op] if op != 'SM' else cache['w_SM']
    mll = cache['mll']
    h, _ = np.histogram(mll, bins=edges, weights=w)
    n, _ = np.histogram(mll, bins=edges)
    return h / h.sum(), n

sp, n_prop = shape_and_n(prop, OP)
sb, n_base = shape_and_n(base, OP)

# relative Poisson error on each shape bin, approximated from raw counts
sig_p = sp / np.sqrt(np.maximum(n_prop, 1))
sig_b = sb / np.sqrt(np.maximum(n_base, 1))
sigma = np.sqrt(sig_p**2 + sig_b**2)

pull = (sp - sb) / sigma
chi2 = np.sum(pull**2)
ndof = len(pull)

print(f"Operator: {OP}")
print(f"{'mll bin':>14} {'pull':>8}")
for i in range(len(edges)-1):
    lo, hi = edges[i], edges[i+1]
    flag = "  <<<" if abs(pull[i]) > 3 else ("  <" if abs(pull[i]) > 2 else "")
    print(f"{lo:6.1f}-{hi:<6.1f} {pull[i]:8.2f}{flag}")

print(f"\nchi2 = {chi2:.1f} over {ndof} bins  ->  chi2/dof = {chi2/ndof:.2f}")
print("(chi2/dof ~ 1 means consistent with pure statistical noise;")
print(" chi2/dof >> 1 means a real shape difference beyond noise)")
