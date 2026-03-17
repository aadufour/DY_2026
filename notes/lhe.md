LHE syntax: drell yan example (p p > z > e+ e-)

#-----------------------------
#INIT BLOCK
<init>
2212 2212 6.500000e+03 6.500000e+03 0 0 247000 247000 -4 1
6.373900e+02 1.338800e+00 6.373900e+02 1
<generator name='MadGraph5_aMC@NLO' version='3.7.0'>please cite 1405.0301 </generator>
</init>

    2212 is the number associated with protons.
    6500 is the energy (in GeV) of the beams: two 6.5 TeV beams = 13 TeV LHC
    0 0 means both beams are unpolarised
    247000  247000 — LHAPDF group ID for the PDF set used (NNPDF2.3 LO) ...pdf stuff
    -4 — weighting strategy. -4 means events are unweighted and the weight is in pb
    1 — number of processes (just one: qq̄ → Z → e⁺e⁻)

6.373900e+02  1.338800e+00  6.373900e+02  1

    This is one entry per process (you'd have more lines if you had more processes):
    637.39 — total cross section in pb
    1.3388 — statistical uncertainty on the cross section (pb)
    637.39 — maximum event weight (same as cross section for unweighted events)
    1 — process ID number

#----------------------
#EVENT
Event header
5 1 +6.3739000e+02 0.7474288E+02 0.7546771E-02 0.1344669E+00

    Event header:
    5 — 5 particles in this event
    1 — process ID 1
    637.39 — event weight (pb)
    74.74 — scale Q (GeV) \sim the Z mass for this event
    0.007547 — α_QED
    0.1345 — α_QCD

Particle 2 — incoming u quark:
2   -1    0    0  501    0   0.0   0.0   137.04   137.04   0.0   0.   1.

    2 — u quark
    -1 — incoming
    0  0 — no mothers (beam particle)
    501  0 — carries color 501, no anticolor
    px=0, py=0 — along beam axis, no transverse momentum
    pz=+137.04 — moving in +z direction (beam 1)
    E=137.04 — massless so E=|pz|
    m=0 — massless
    0. — lifetime
    1. — helicity +1/2

    x_1 = 137.04 / 6500 = 0.02108 ← you'll see this again in <pdfrwt beam="1">

Particle 2 — incoming ubar antiquark:
-2   -1    0    0    0  501   0.0   0.0   -10.19   10.19   0.0   0.  -1.

    -2 — ubar
    -1 — incoming
    0  501 — carries anticolor 501, no color (color-matched with the u quark above)
    pz=-10.19 — moving in -z direction (beam 2)
    E=10.19 — again E=|pz|
    -1. — helicity -1/2

    x₂ = 10.19 / 6500 = 0.001568 ← again matches <pdfrwt beam="2">
    Notice x_1 >> x_2, so the Z will be strongly boosted in the +z direction.

Particle 3 — Z boson:
23    2    1    2    0    0   0.0   0.0   126.85   147.23   74.74   0.   9.

    23 — Z boson
    2 — intermediate (status=2, it will decay)
    1  2 — mothers are particles 1 and 2 (the u and ubar) ✓
    0  0 — colorless
    px=0, py=0 — no transverse momentum (LO: initial quarks have none)
    pz=+126.85 — strongly boosted in +z, as expected from x_1>>x_2
    E=147.23 — total energy
    m=74.74 — off-shell! Z pole is 91.19 GeV but this one is at 74.74 GeV, sampled from the Breit-Wigner
    9. — helicity label 9 = undefined/mixed for a decaying boson

    You can verify 4-momentum conservation:

    pz: 137.04 + (-10.19) = 126.85 
    E: 137.04 + 10.19 = 147.23

Particle 4 — outgoing e+ (positron):
-11    1    3    3   0   0   20.79  -14.75   9.60   27.24   0.0   0.  -1.

    -11 — positron
    1 — final state
    3  3 — mother is particle 3 (the Z) 
    px=+20.79, py=-14.75, pz=+9.60
    E=27.24
    m=0 — massless
    -1. — helicity -1/2

Particle 5 — outgoing e- (electron):
11    1    3    3   0   0  -20.79  +14.75   117.25   120.00   0.0   0.   1.

    11 — electron
    1 — final state
    3  3 — also from the Z 
    px=-20.79, py=+14.75 — exactly opposite to the positron in x and y  (momentum conservation)
    pz=+117.25 — most of the Z's boost goes to the electron here
    E=120.00
    1. — helicity +1/2

    Check lepton momentum conservation against the Z:

    px: 20.79 + (-20.79) = 0 
    py: -14.75 + 14.75 = 0 
    pz: 9.60 + 117.25 = 126.85 
    E: 27.24 + 120.00 = 147.24

#<mgrwt> BLOCK
<rscale>  0  74.74 </rscale>

    0 — zero QCD vertices (LO, no α_s in this diagram)
    74.74 — renormalization scale used (= Z mass for this event)

<asrwt>0</asrwt>

    alpha_s reweighting power = 0, consistent with LO electroweak process

<pdfrwt beam="1">  1  2  0.02108  74.74 </pdfrwt>
<pdfrwt beam="2">  1 -2  0.001568 74.74 </pdfrwt>

    1 — PDF set index
    2 / -2 — parton flavor (u / ubar)
    0.02108 — x₁ (= 137.04/6500) 
    0.001568 — x₂ (= 10.19/6500) 
    74.74 — factorization scale Q at which PDFs were evaluated

<totfact> 23128.9 </totfact>
    This is the combined factor: PDF₁(x₁,Q) × PDF₂(x₂,Q) × matrix element / phase space. 
    It's used internally by MadGraph to reweight events for systematics variations (different scales or PDF sets) without regenerating them.