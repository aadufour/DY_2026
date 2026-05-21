#!/bin/bash
# LLR morphing environment setup.
# Uses CMSSW_spritz with AnalyticAnomalousCoupling (template_morphing) + tools (morphing_model)

CMSSW_DIR="/grid_mnt/data__data.polcms/cms/adufour/CMSSW_spritz/CMSSW_14_1_0_pre4/src"
TOOLS_DIR="${CMSSW_DIR}/tools"

cd "${CMSSW_DIR}" && cmsenv && cd - > /dev/null
cd "${TOOLS_DIR}" && source env.sh && cd - > /dev/null

# Prepend morphing tools to PATH so they take precedence over repo's combine_tools.
# Order matters: combine_helpers (createWS, runScans) > combination (createCombineJson, createJson) > plotters
export PATH="${TOOLS_DIR}/combine_helpers:${TOOLS_DIR}/combination:${TOOLS_DIR}/plotters:${PATH}"

echo "LLR combine MORPHING environment ready."
echo "  CMSSW : ${CMSSW_DIR}"
echo "  tools : ${TOOLS_DIR} (morphing_model)"
echo "  PATH  : combine_helpers + combination + plotters prepended"
