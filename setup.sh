#therelease=CMSSW_9_4_9
#export SCRAM_ARCH=slc6_amd64_gcc630
#mkdir -p common
#if [ ! -d common/$therelease ]; then 
#    cd common/ ;
#    cmsrel $therelease;
#    cd ~-
#fi
#cd common/$therelease/src
#eval `scram ru -sh`
#cd -
#

export FTBASE="$( cd "$(dirname "$BASH_SOURCE")" ; pwd -P )"

echo $FTBASE
#export LD_LIBRARY_PATH=${FTBASE}/babymaking/batch/:$LD_LIBRARY_PATH
export PYTHONPATH=${FTBASE}/:$PYTHONPATH
export PYTHONPATH=${FTBASE}/common/matplottery:$PYTHONPATH

# bash PS1 gets mangled because CMSSW overrides this to LANG=C :(, so switch it back.
export LANG=en_US.UTF-8

# export PYTHONPATH=$PWD/analysis/bdt/xgboost/python-package/lib/python2.7/site-packages/:$PYTHONPATH
# export PYTHONPATH=$PWD/analysis/bdt/xgboost/python-package/:$PYTHONPATH
# export PYTHONPATH=$PWD/analysis/bdt/root_numpy-4.7.2/lib/python2.7/site-packages/:$PYTHONPATH

[[ -d ${FTBASE}/common/matplottery/ ]] || {
    git clone git@github.com:aminnj/matplottery.git ${FTBASE}/common/matplottery/;
    pip install --user matplotlib
    pip install --user uproot
    export PATH=/home/users/${USER}/.local/bin:$PATH
}
# [[ -d babymaking/batch/NtupleTools/ ]] || git clone git@github.com:cmstas/NtupleTools.git babymaking/batch/NtupleTools/
#if [ ! -d ${FTBASE}/babymaking/batch/ProjectMetis/ ]; then
#    git clone https://github.com/aminnj/ProjectMetis/ ${FTBASE}/babymaking/batch/ProjectMetis/
#else
#    source ${FTBASE}/babymaking/batch/ProjectMetis/setup.sh
#fi
## [[ -d ${FTBASE}/common/Software/ ]] || git clone git@github.com:cmstas/Software.git ${FTBASE}/common/Software/
[[ -d ${FTBASE}/common/CORE/ ]] || {
    git clone git@github.com:cmstas/CORE.git ${FTBASE}/common/CORE/;
    cd ${FTBASE}/common/CORE; 
    make -j10 >&/dev/null &
    cd -;
}

#[[ -d ${FTBASE}/common/${therelease}/src/HiggsAnalysis/CombinedLimit/ ]] || {
#    pushd ${FTBASE}/
#    cd $CMSSW_BASE/src
#    cmsenv
#    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
#    cd HiggsAnalysis/CombinedLimit
#   git checkout 94x
#    scramv1 b vclean
#    scramv1 b -j 15
#    popd
#}
#
