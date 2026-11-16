source ~/.bashrc; dy_combine_morphing
OPT="-t -1 -m 125 -v 0 --robustFit=1 --setRobustFitTolerance=0.2 --cminDefaultMinimizerStrategy=0 --X-rtd=MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=99999999999 --cminFallbackAlgo Minuit2,Migrad,0:0.2 --stepSize=0.005 --X-rtd FITTER_NEW_CROSSING_ALGO --X-rtd FITTER_NEVER_GIVE_UP --X-rtd FITTER_BOUND"
C=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs
for v in propcorr_v1_minnlo propcorr_v3_pdfcut; do
  d=$C/$v/datacards/inc_mm/mll; mkdir -p $d/pdf_constraint; cd $d/pdf_constraint
  for op in cHQ1 cHl1; do
    r=$(python3 -c "import json; a,b=json.load(open(\"../metadata.json\"))[\"operators\"][\"$op\"]; print(f\"{a},{b}\")")
    combine -M MultiDimFit ../model_$op.root --algo singles --redefineSignalPOIs PDFweight -P PDFweight --setParameters r=1,k_$op=0 --freezeParameters r,k_$op --setParameterRanges PDFweight=-3,3 $OPT -n .pdfsingles_$op > singles_$op.log 2>&1 &
    combine -M MultiDimFit ../model_$op.root --algo grid --points 100 --redefineSignalPOIs k_$op -P k_$op --setParameters r=1,k_$op=0,PDFweight=0 --freezeParameters r,PDFweight --setParameterRanges k_$op=$r $OPT -n .pdffrozen_$op > frozen_$op.log 2>&1 &
  done
done
wait; echo CONSTRDONE
