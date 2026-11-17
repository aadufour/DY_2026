#!/usr/bin/env python3
"""
How much of each operator's linear effect is just a normalisation change, and how degenerate
pairs become once the normalisation is free (as with lumi / cross-section nuisances).

Per variable, with stat weights 1/sqrt(sm+bkg) per bin:
  |cos(lin_op, sm)|  : 1 = the operator only rescales the SM shape (degenerate with a norm. nuisance)
  sqrt(1-cos^2)      : fraction of the sensitivity left once the normalisation floats
  axis ratio         : elongation of the linear 2D ellipse after projecting out the SM direction,
                       sqrt((1+|c|)/(1-|c|)) with c the cosine of the projected lin templates
Run from the configs dir (reads propcorr_v3_2d/datacards/inc_mm/<var>/shapes.root).
"prob" = pairs flagged by hand as bad 2D mll scans (2026-11-17).
"""
import ROOT, numpy as np, itertools
prob={frozenset(p.split("_")) for p in "cHQ1_cHbq cHQ1_cHd cHQ1_cHu cHQ1_cll1 cHbq_cHl3 cHbq_cll1 cHd_cHbq cHd_cHl3 cHj1_cHQ1 cHj1_cHbq cHj1_cHd cHj1_cHe cHj1_cHl1 cHj1_cHu cHj1_cQe cHj3_cHbq cHj3_cHu cHu_cHbq cHu_cHd cHu_cll1 cbBRe_cHbq cbBRe_cHd cbBRe_cHj1 cbWRe_cHbq cbWRe_cHd cbWRe_cHj1 cbWRe_cbBRe".split()}
for var in ["mll","triple_diff"]:
    f=ROOT.TFile.Open(f"propcorr_v3_2d/datacards/inc_mm/{var}/shapes.root"); names={k.GetName() for k in f.GetListOfKeys()}
    A=lambda n: np.array([f.Get(n).GetBinContent(i) for i in range(1,f.Get(n).GetNbinsX()+1)])
    sm=A("histo_sm"); bkg=sum(A(n) for n in names if n.startswith("histo_") and not any(t in n for t in ("Up","Down","Data","sm","w1_","wm1_","w11_")))
    sig=np.sqrt(sm+bkg); s=sm/sig; s/=np.linalg.norm(s)
    ops=[o for o in sorted(n[9:] for n in names if n.startswith("histo_w1_") and not n.endswith(("Up","Down"))) if o not in ("cHQ3","cQl3")]
    L={o:(A(f"histo_w1_{o}")-A(f"histo_wm1_{o}"))/2/sig for o in ops}
    print(f"== {var}: |cos(lin, SM shape)| = fraction of the linear effect that is a pure normalisation change")
    cs={o:abs(L[o]@s)/np.linalg.norm(L[o]) for o in ops}
    print("   "+"  ".join(f"{o}:{cs[o]:.3f}" for o in sorted(ops,key=lambda o:-cs[o])))
    print(f"   residual sensitivity after floating the normalisation sqrt(1-cos^2): " + "  ".join(f"{o}:{np.sqrt(1-cs[o]**2):.2f}" for o in sorted(ops,key=lambda o:-cs[o])[:10]))
    R={o:L[o]-(L[o]@s)*s for o in ops}
    rat={}
    for a,b in itertools.combinations(ops,2):
        c=R[a]@R[b]/np.linalg.norm(R[a])/np.linalg.norm(R[b]); rat[frozenset((a,b))]=np.sqrt((1+abs(c))/max(1-abs(c),1e-12))
    P=[rat[k] for k in rat if k in prob]; O=[rat[k] for k in rat if k not in prob]
    print(f"   axis ratio with free normalisation: problem median {np.median(P):.1f} [{np.percentile(P,10):.1f},{np.percentile(P,90):.1f}]   others median {np.median(O):.1f} [{np.percentile(O,10):.1f},{np.percentile(O,90):.1f}]")
    hit=sorted(rat,key=lambda k:-rat[k])[:40]; print(f"   of the 40 most elongated pairs, {sum(k in prob for k in hit)} are in the problem list (27 of {len(rat)} pairs overall)")
