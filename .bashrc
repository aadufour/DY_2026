
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/usr/local/envs/spritz/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/local/envs/spritz/etc/profile.d/conda.sh" ]; then
        . "/usr/local/envs/spritz/etc/profile.d/conda.sh"
    else
        export PATH="/usr/local/envs/spritz/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

source /cvmfs/cms.cern.ch/cmsset_default.sh
export MG5_PATH=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/madgraph5amcatnlo/2.7.3-3906416aae233758192cea04e92e2d24
alias mg5="$MG5_PATH/bin/mg5_aMC"
export LD_LIBRARY_PATH=/home/llr/cms/adufour/MG5/mg5amcnlo/HEPTools/lhapdf6_py3/lib:$LD_LIBRARY_PATH
export PATH=/home/llr/cms/adufour/mg5-fortran12/bin:$PATH
export PATH=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28/bin:$PATH
export LD_LIBRARY_PATH=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/gcc/12.3.1-40d504be6370b5a30e3947a6e575ca28/lib:$LD_LIBRARY_PATH
alias cdwork='cd /grid_mnt/data__data.polcms/cms/adufour'


export PATH="$PATH:/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools"

dy_analysis() {
    source /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools/env_llr.sh analysis
}
dy_combine() {
    source /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools/env_llr.sh
}
dy_combine_morphing() {
    source /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools/env_llr_morphing.sh
}


export EOS_MGM_URL=root://eosuser.cern.ch

alias t3stat='/opt/exp_soft/cms/t3/t3stat'

export PATH=$PATH:/grid_mnt/data__data.polcms/cms/adufour/CMSSW_14_1_0_pre4/src/tools/systematic_tools
export PATH=/grid_mnt/data__data.polcms/cms/adufour/spritz/analysis/spritz:$PATH
export PATH=/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/spritz:$PATH
export PATH=/grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools:$PATH
export PYTHONPATH=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/src:$PYTHONPATH
export SPRITZ_PATH=/grid_mnt/data__data.polcms/cms/adufour/spritz_fabian

runScans.py() { python3 /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools/runScans.py "$@"; }
runPlots.py() { python3 /grid_mnt/data__data.polcms/cms/adufour/DY_2026/analysis/combine_tools/runPlots.py "$@"; }
export -f runScans.py runPlots.py



alias spritz-shell='apptainer exec -B /etc/grid-security/certificates:/etc/grid-security/certificates -B /cvmfs -B /grid_mnt -B /grid_mnt/data__data.polcms/cms/adufour/spritz_fabian/data/Full2018v9/samples/samples.json:/opt/spritz/data/Full2018v9/samples/samples.json /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif bash --rcfile ~/.bashrc'


#tmux stuff

tlog() {
    local name=${1:-work}
    local log="$PWD/tmux-${name}.log"
    echo "===== session '$name' started $(date '+%Y-%m-%d %H:%M:%S') =====" >> "$log"
    tmux new-session -A -s "$name" \; pipe-pane -o "cat >> $log"
}

alias spritz-shell-giacomo='PYTHONPATH=/grid_mnt/data__data.polcms/cms/adufour/spritz_giacomo/src:$PYTHONPATH apptainer exec -B /etc/grid-security/certificates:/etc/grid-security/certificates -B /cvmfs -B /grid_mnt -B /grid_mnt/data__data.polcms/cms/adufour/spritz_giacomo/data/Full2018v9/samples/samples.json:/opt/spritz/data/Full2018v9/samples/samples.json /grid_mnt/data__data.polcms/cms/adufour/spritz-env.sif bash --rcfile ~/.bashrc'
