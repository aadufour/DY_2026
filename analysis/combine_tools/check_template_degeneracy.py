# Find degenerate operator pairs from the datacard templates (w1/wm1/sm): bit-identical lin+quad, or |cos(lin)|>0.99.
# Run from datacards/inc_mm (needs mll/metadata.json and <var>/shapes.root).
import ROOT, numpy as np, json, itertools, sys
ops = list(json.load(open("mll/metadata.json"))["operators"])
for var in ["triple_diff", "mll"]:
    f = ROOT.TFile.Open(f"{var}/shapes.root")
    names = [k.GetName() for k in f.GetListOfKeys()]
    if var == "triple_diff":
        print([n for n in names if "Up" not in n and "Down" not in n and "w11" not in n][:12])
    def arr(n):
        h = f.Get(n); return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)])
    pre = [n for n in names if n.endswith("_sm") and "Up" not in n and "Down" not in n][0][:-2]
    sm = arr(pre + "sm")
    lin, quad = {}, {}
    for o in ops:
        wp, wm = arr(pre + "w1_" + o), arr(pre + "wm1_" + o)
        lin[o] = (wp - wm) / 2; quad[o] = (wp + wm) / 2 - sm
    print(f"== {var}: identical (lin & quad bit-equal) / cos(lin) > 0.99 ==")
    for a, b in itertools.combinations(ops, 2):
        same = np.array_equal(lin[a], lin[b]) and np.array_equal(quad[a], quad[b])
        na, nb = np.linalg.norm(lin[a]), np.linalg.norm(lin[b])
        c = lin[a] @ lin[b] / (na * nb) if na > 0 and nb > 0 else 0
        cq = quad[a] @ quad[b] / (np.linalg.norm(quad[a]) * np.linalg.norm(quad[b]))
        if same or abs(c) > 0.99:
            print(f"  {a:6s} {b:6s} identical={same} cos_lin={c:+.5f} cos_quad={cq:+.5f} |lin| ratio={na/nb:.4f}")
