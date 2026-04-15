#!/bin/bash
# Source this before running any combine scripts on LLR:
#   cd /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine
#   source env_llr.sh

# 1. CMSSW environment (provides HiggsAnalysis.CombinedLimit.*)
CMSSW_DIR="/grid_mnt/data__data.polcms/cms/adufour/CMSSW_14_1_0_pre4/src"
cd "${CMSSW_DIR}" && cmsenv && cd - > /dev/null

# 2. combine_tools (provides python.DatacardHelpers, plotters, combine_helpers, etc.)
TOOLS_DIR="/grid_mnt/data__data.polcms/cms/adufour/CMSSW_14_1_0_pre4/src/tools"
cd "${TOOLS_DIR}" && source env.sh && cd - > /dev/null

echo "LLR combine environment ready."
echo "  CMSSW : ${CMSSW_DIR}"
echo "  tools : ${TOOLS_DIR}"
