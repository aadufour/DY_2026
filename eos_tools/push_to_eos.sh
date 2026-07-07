#!/bin/bash
set -e

CONFIG=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/configs/eft_bkg_fullsyst_v6
INC_MM=$CONFIG/datacards/inc_mm
PLOTS=$CONFIG/plots

EOS_XROOTD_BASE=root://eosuser.cern.ch//eos/user/a/aldufour/www/
EOS_SUBDIR=RECO3D_v2

# adjust this if pushtoEOS_parallel.py doesn't live in your $PWD/$PATH
PUSH_SCRIPT=pushtoEOS_parallel.py


for sub in mll rapll_abs costhetastar triple_diff; do
    echo "=== Processing $sub ==="

    for kind in scan nuis_eft nuis_combine; do
        echo "--> uploading $sub/$kind"
        python3 "$PUSH_SCRIPT" "$INC_MM/$sub/$kind" "$EOS_SUBDIR/$sub/$kind" -j 32
    done

    echo "--> uploading $sub dist (from plots/$sub)"
    python3 "$PUSH_SCRIPT" "$PLOTS/$sub" "$EOS_SUBDIR/$sub/dist" -j 32

    echo "--> uploading $sub eft_summary_two_panel.pdf/png"
    xrdcp -f "$INC_MM/$sub/eft_summary_two_panel.pdf" "$EOS_XROOTD_BASE/$EOS_SUBDIR/$sub/"
    xrdcp -f "$INC_MM/$sub/eft_summary_two_panel.png" "$EOS_XROOTD_BASE/$EOS_SUBDIR/$sub/"


done


echo "--> uploading top-level eft_summary_multivar_horizontal.pdf/png"
xrdcp -f "$CONFIG/eft_summary_multivar_horizontal.pdf" "$EOS_XROOTD_BASE/$EOS_SUBDIR/"
xrdcp -f "$CONFIG/eft_summary_multivar_horizontal.png" "$EOS_XROOTD_BASE/$EOS_SUBDIR/"


echo "All done."
