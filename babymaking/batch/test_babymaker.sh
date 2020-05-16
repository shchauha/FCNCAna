#!/usr/bin/env bash

function get_dataset() {
    which=$1
    sampstr=$2
    python -c 'from samples import '${which}' ; y = [x for x in '${which}' if "'${sampstr}'" in x[0]]; print y[len(y)//2][0]';
}

function get_random_file() {
    datasetname=$1;
    dir=$(dis_client.py -t snt "${datasetname} | grep location" | head -1);
    last=$(ls -1 ${dir}/*.root | sort | shuf --random-source=<(yes) | head -n 1)
    echo $last
}

function get_file() {
    get_random_file $(get_dataset $1 "$2")
}

function check() {
    fname=$1
    badstr="\e[91mx\e[0m"
    goodstr="\e[92mo\e[0m"
    pgood=$([ -n "$(grep "Processed [1-9]" ${fname})" ] && echo $goodstr || echo $badstr);
    sgood=$([ -n "$(grep "Stored [1-9]" ${fname})" ] && echo $goodstr || echo $badstr);
    echo -e "$pgood$sgood $fname"
}

function run() {
    echo Running ./main.exe $(get_file $1 $2) -o $3 -n 50000 '>&' $4
    ./main.exe $(get_file $1 $2) -o $3 -n 50000 >& $4
    check $4
}

outputdir=outputs_test/
mkdir -p ${outputdir}
rm ${outputdir}/log*.txt

echo $(get_random_file $(get_dataset mc_2017 "/TTWJets"))
return

# run mc_2016       "/TTWJets"    ${outputdir}/test_2016_ttw.root        ${outputdir}/log_2016_ttw.txt        &
# run mc_2016       "/TTTT"       ${outputdir}/test_2016_tttt.root       ${outputdir}/log_2016_tttt.txt       &
run mc_2016_94x       "T1tttt"       ${outputdir}/test_2016_94x_smst1tttt.root       ${outputdir}/log_2016_94x_smst1tttt.txt       &
run mc_2016_94x   "/TTTT"       ${outputdir}/test_2016_94x_tttt.root   ${outputdir}/log_2016_94x_tttt.txt   &
# run data_2016     "/DoubleEG"   ${outputdir}/test_2016_dataee.root     ${outputdir}/log_2016_dataee.txt     &
# run data_2016     "/DoubleMuon" ${outputdir}/test_2016_datamm.root     ${outputdir}/log_2016_datamm.txt     &
run data_2016_94x "/DoubleEG"   ${outputdir}/test_2016_94x_dataee.root ${outputdir}/log_2016_94x_dataee.txt &
run data_2016_94x "/DoubleMuon" ${outputdir}/test_2016_94x_datamm.root ${outputdir}/log_2016_94x_datamm.txt &
run mc_2017       "/TTWJets"    ${outputdir}/test_2017_ttw.root        ${outputdir}/log_2017_ttw.txt        &
run mc_2017       "T1tttt"       ${outputdir}/test_2017_smst1tttt.root       ${outputdir}/log_2017_smst1tttt.txt       &
# run mc_2017       "/ttHToNonbb" ${outputdir}/test_2017_tth.root        ${outputdir}/log_2017_tth.txt        &
run mc_2017       "/TTTT"       ${outputdir}/test_2017_tttt.root       ${outputdir}/log_2017_tttt.txt       &
run data_2017     "/DoubleEG"   ${outputdir}/test_2017_dataee.root     ${outputdir}/log_2017_dataee.txt     &
run data_2017     "/DoubleMuon" ${outputdir}/test_2017_datamm.root     ${outputdir}/log_2017_datamm.txt     &
run data_2018     "/EGamma"     ${outputdir}/test_2018_dataee.root     ${outputdir}/log_2018_dataee.txt     &
run data_2018     "/DoubleMuon" ${outputdir}/test_2018_datamm.root     ${outputdir}/log_2018_datamm.txt     &
run mc_2018       "/TTTT"       ${outputdir}/test_2018_tttt.root       ${outputdir}/log_2018_tttt.txt       &
run mc_2018       "T1tttt"       ${outputdir}/test_2018_smst1tttt.root       ${outputdir}/log_2018_smst1tttt.txt       &
# run data_2018     "17Sep2018-v1" ${outputdir}/test_2018_datarereco.root     ${outputdir}/log_2018_datarereco.txt     &

wait
