#!/bin/bash
# LLR morphing environment setup.
# Uses CMSSW_spritz with AnalyticAnomalousCoupling (template_morphing) + tools (morphing_model)

CMSSW_DIR="/grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4/src"
TOOLS_DIR="${CMSSW_DIR}/tools"

cd "${CMSSW_DIR}" && cmsenv && cd - > /dev/null
cd "${TOOLS_DIR}" && source env.sh && cd - > /dev/null

# Morphing tools on PATH (runPlots_compare, mkEFTScan, ...), with the repo's combine_tools FIRST:
# runScans.py, createWS.py, createCombineJson.py, createJson.py, runPlots.py exist in both, and the
# repo versions are the maintained ones (the tools copies date from 2026-05 and lack e.g. --points,
# --condor, --ops). Shell functions in ~/.bashrc already pointed to the repo versions, but nice/xargs/
# condor/subprocess bypass functions and used to pick the stale tools copies.
REPO_TOOLS="/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools"
export PATH="${REPO_TOOLS}:${TOOLS_DIR}/combine_helpers:${TOOLS_DIR}/combination:${TOOLS_DIR}/plotters:${PATH}"

echo "LLR combine MORPHING environment ready."
echo "  CMSSW : ${CMSSW_DIR}"
echo "  tools : ${TOOLS_DIR} (morphing_model)"
echo "  PATH  : repo combine_tools, then combine_helpers + combination + plotters, prepended"
