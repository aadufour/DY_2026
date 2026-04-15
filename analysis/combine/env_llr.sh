#!/bin/bash
# LLR environment setup.
#
# Two modes:
#   source env_llr.sh          → CMSSW (for createWS, runScans, createCombineJson)
#   source env_llr.sh analysis → analysis venv (for build_cache, build_datacard)
#
# Workflow:
#   Step 1  build_cache.py       → source env_llr.sh analysis
#   Step 2  build_datacard.py    → source env_llr.sh analysis
#   Step 3  createJson.py        → either
#   Step 4  createCombineJson.py → source env_llr.sh
#   Step 5  createWS.py          → source env_llr.sh
#   Step 6  runScans.py          → source env_llr.sh
#   Step 7  runPlots.py          → source env_llr.sh

ANALYSIS_VENV="/grid_mnt/data__data.polcms/cms/adufour/analysis_venv"
CMSSW_DIR="/grid_mnt/data__data.polcms/cms/adufour/CMSSW_14_1_0_pre4/src"
TOOLS_DIR="${CMSSW_DIR}/tools"

if [[ "$1" == "analysis" ]]; then
    source "${ANALYSIS_VENV}/bin/activate"
    echo "Analysis environment ready (numpy $(python3 -c 'import numpy; print(numpy.__version__)'))."
    echo "  venv : ${ANALYSIS_VENV}"
else
    cd "${CMSSW_DIR}" && cmsenv && cd - > /dev/null
    cd "${TOOLS_DIR}" && source env.sh && cd - > /dev/null
    echo "LLR combine environment ready."
    echo "  CMSSW : ${CMSSW_DIR}"
    echo "  tools : ${TOOLS_DIR}"
fi
