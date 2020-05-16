import sys
import glob
import os
import time
import itertools
import json
import traceback

from metis.Sample import SNTSample,DirectorySample
from metis.LocalMergeTask import LocalMergeTask
from metis.CondorTask import CondorTask
from metis.StatsParser import StatsParser
from metis.Utils import send_email

def get_event_rate(fname):
    fname = fname.replace(".err", ".out")
    rate = -1
    if not os.path.exists(fname): return rate
    with open(fname, "r") as fhin:
        for line in fhin:
            if line.startswith("Processed "):
                parts = line.split()
                den = float(parts[4])
                if den < 1e-6: break
                rate = 1.0*int(parts[1])/den
                break
    return rate

def split_func(dsname):
    if any(x in dsname for x in [
        "/W","/Z","/TTJets","/DY","/ST",
        ]):
        return 5
    elif "/EGamma" in dsname:
        return 12
    elif "Run201" in dsname:
        return 7
    elif "/SMS-" in dsname:
        return 2
    elif "LOoblique" in dsname:
        return 30
    else:
        return 3

if __name__ == "__main__":

    from samples import data_2016, mc_2016, data_2017, mc_2017, data_2018, mc_2018, data_2016_94x, mc_2016_94x

    # year_sample_map = [("2016",data_2016+mc_2016)] + [("2017",data_2017+mc_2017)] + [("2018",data_2018)] + [("2016_94x",data_2016_94x+mc_2016_94x)]
    year_sample_map = [("2016",data_2016+mc_2016)] + [("2017",data_2017+mc_2017)] + [("2018",data_2018+mc_2018)]
    # year_sample_map = [("2016_94x",data_2016_94x+mc_2016_94x)]
    # tag = "v3.16_all" # back to old V6 JECs
    # tag = "v3.17_all" # V6 JECs, pdf weight clipping, fixed JetID for jec/jer up/down
    # tag = "v3.19_all" # same as v17 with top tagging vars and 2 more potential bdt variables
    # tag = "v3.19_newdata" # to validate event counts with new prompt 2018Dv2
    # tag = "v3.20_jec32" # try JECv32 and new WPs with 2017, 2018 only
    # tag = "v3.21" # v32 JECs. all years, new 10-02-04 MC, including some 2018 MC -- bugged 2017 missing 2 ele iso branches
    # year_sample_map =  [("2017",data_2017+mc_2017)] # NOTE
    # tag = "v3.21_fix2017" # 2017 only. fixing iso issues
    # year_sample_map =  [("2016",data_2016+mc_2016)]+[("2018",data_2018+mc_2018)] # NOTE
    # year_sample_map =  [("2018",data_2018+mc_2018)] # NOTE
    # year_sample_map =  [("2016",data_2016+mc_2016)] # NOTE
    # year_sample_map =  [("2016",mc_2016)] # NOTE
    # year_sample_map = [("2017",data_2017+mc_2017)]
    # year_sample_map = [("2018",data_2018+mc_2018)]
    # tag = "v3.22" # 2018 only. latest ID (v2), remade corruptions, all MC samples included I think
    # tag = "v3.23" # 2016,2018 only. for 2018, latest ID (v2), remade corruptions, all MC samples included I think -- latest 102X btag WPs
    # tag = "v3.23_oldbwp2018" # v3.23 with old (matching 2017) btag WP for 2018

    #year_sample_map = [
    #        ("2016",mc_2016),
    #        ("2017",mc_2017),
    #        ]
    ##,("2018",mc_2018)]
    ## year_sample_map = [("2016",mc_2016)]
    ## year_sample_map = [("2017",mc_2017)]
    ## year_sample_map = [("2016",mc_2016),("2017",mc_2017)]
    #tag = "v3.24" # all 3 years. split btag unc, fastsim 2016/2017
    #extra_args = ""
    #tag_match = ""
    ## remember to replace in samples after it finishes
    ## year_sample_map = [("2016",[
    ##     ["/SMS-T6ttHZ_BR-H_0p6_mStop300to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM|/hadoop/cms/store/group/snt/run2_mc2016_94x/SMS-T6ttHZ_BR-H_0p6_mStop300to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_MINIAODSIM_CMS4_V10-02-05/", "T6TTHZ_m1"],
    ##     ["/SMS-T6ttHZ_BR-H_0p6_mStop1050to1600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM|/hadoop/cms/store/group/snt/run2_mc2016_94x/SMS-T6ttHZ_BR-H_0p6_mStop1050to1600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv3-PUSummer16v3Fast_94X_mcRun2_asymptotic_v3-v1_MINIAODSIM_CMS4_V10-02-05/", "T6TTHZ_m2"],
    ##     ])]
    ## extra_ags = "--ignorebadfiles"

    # year_sample_map = [("2018",mc_2018)]
    # # year_sample_map = [("2018",data_2018)]
    # # year_sample_map = [("2018",data_2018)]
    # # tag = "v3.25" # 2018 only, preliminary jecs without residuals, before they were announced (so they could be wrong) -- bugged Autumn18_V1 JECs (but identical to V2 actually)
    # # tag = "v3.26" # 2018 only, preliminary jecs without residuals -- bugfixed Autumn18_V3 JECs
    # # tag = "v3.27" # 2018 only, same as v3.26 with new 102x btag SFs
    # # tag = "v3.27_resid17" # same as v3.27 (but DATA ONLY, so need to copy MC over) but fall17F V32 residuals (Autumn18_V3_OldRes data payload)
    # tag = "v3.28" # new V8 resids split into ABCD
    # extra_args = ""
    # tag_match = ""
    # # extra_args = "--xrootd" # cabinet nodes have messed up hadoop mapping
    # # extra_args = "--ignorebadfiles" # FIXME

    # year_sample_map = [
    #         ("2017",mc_2017),
    #         ]
    # tag = "v3.24_newfixfsjecs" # 2017 fastsim jecs
    # tag_match = ""
    # year_sample_map = [("2017",[
    #     ["/SMS-T1tttt_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_pilot_94X_mc2017_realistic_v15-v3/MINIAODSIM", "T1TTTT"],
    #     ["/SMS-T1tttt_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_pilot_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T1TTTT_ext"],
    #     ["/SMS-T6ttWW_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T6TTWW"],
    #     ["/SMS-T6ttWW_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T6TTWW_ext"],
    #     ["/SMS-T5qqqqVV_dM20_mGlu-600to2300_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T5QQQQVV_dm20"],
    #     ["/SMS-T5qqqqVV_dM20_mGlu-600to2300_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T5QQQQVV_dm20_ext"],
    #     ["/SMS-T5qqqqVV_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T5QQQQVV_main"],
    #     ["/SMS-T5qqqqVV_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T5QQQQVV_main_ext"],
    #     ["/SMS-T1ttbb_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T1TTBB"],
    #     ["/SMS-T5tttt_dM175_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T5TTTT"],
    #     ["/SMS-T5tttt_dM175_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T5TTTT_ext"],
    #     ["/SMS-T1qqqqL_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T1QQQQL_main"],
    #     ["/SMS-T1qqqqL_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T1QQQQL_main_ext"],
    #     ["/SMS-T6ttHZ_BR-H_0p6_mStop300to1000_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T6TTHZ_m1"],
    #     ["/SMS-T6ttHZ_BR-H_0p6_mStop1050to1600_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15-v1/MINIAODSIM", "T6TTHZ_m2"], # Doesn't overlap with 2016 sample
    #     ["/SMS-T6ttHZ_BR-H_0p6_mStop300to1000_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T6TTHZ_m1_ext"],
    #     ["/SMS-T6ttHZ_BR-H_0p6_mStop1050to1600_TuneCP2_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PUFall17Fast_94X_mc2017_realistic_v15_ext1-v1/MINIAODSIM", "T6TTHZ_m2_ext"], # Doesn't overlap with 2016 sample
    #     ])]
    # extra_args = "--ignorebadfiles"

    # year_sample_map = [("2017",[
    #     ["/WZTo3LNu_0Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "test_WZmll4_0j"],
    #     ["/WZTo3LNu_1Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "test_WZmll4_1j"],
    #     ["/WZTo3LNu_2Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "test_WZmll4_2j"],
    #     ["/WZTo3LNu_3Jets_MLL-4to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "test_WZmll4_3j"],
    #     ["/WZTo3LNu_0Jets_MLL-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/MINIAODSIM", "test_WZmll50_0j"],
    #     ["/WZTo3LNu_1Jets_MLL-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM", "test_WZmll50_1j"],
    #     ["/WZTo3LNu_2Jets_MLL-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM", "test_WZmll50_2j"],
    #     ["/WZTo3LNu_3Jets_MLL-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM", "test_WZmll50_3j"],
    #     ])]
    # extra_ags = "--skipos --ignorebadfiles"

    # extra_args = "--ignorebadfiles" # FIXME
    # year_sample_map = [("2016",[
    #     ["/ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", "ST2"],
    #     ["/ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", "ST1"],
    #     ])]


    # extra_args = ""
    # tag = "v3.29" # 94x 2016
    # tag_match = ""
    # year_sample_map = [("2016_94x",data_2016_94x+mc_2016_94x)]

    extra_args = ""
    # tag = "v3.30" # 2018 recovered EGammaDv2
    # tag = "v3.30_hhat"
    tag = "v3.31"
    tag_match = ""
    year_sample_map = [("2016_94x",data_2016_94x+mc_2016_94x)] + [("2017",data_2017+mc_2017)] + [("2018",data_2018+mc_2018)]
    # year_sample_map = [("2018",[
    #     ["/EGamma/Run2018D-22Jan2019-v2/MINIAOD", "DataDoubleEGDv2"],
    #     ])]
    extra_args = ""
    tag = "v3.31_jec15"
    tag_match = ""
    year_sample_map = [("2018",data_2018)]

    # extra_args = ""
    # tag = "v3.31_testc"
    # tag_match = ""
    # year_sample_map = [("2018",[
    #     ["/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM", "TTHtoNonBB"],
    #     ["/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM", "TTWnlo"],
    #     ["/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM", "TTZnlo"],
    #     ])]

    extra_args = "--skipos --ignorebadfiles"
    tag = "v3.31_ttslv2"
    tag_match = ""
    year_sample_map = [
            # ("2016",[
            #     ["/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "TTSLtop"],
            #     ["/TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/MINIAODSIM", "TTSLtopbar"],
            #     ]),
            # ("2017",[
            #     ["/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "TTSLtop"],
            #     ["/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM", "TTSLtopbar"],
            #     ]),
            ("2018",[
                ["/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "TTSLtop"],
                ["/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM", "TTSLtopbar"],
                ["/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM", "TTinclusive"],
                ]),
            ]



    # Loop 'til we get dizzy
    for i in range(1000):
        ndone = 0
        total_summary = {}
        nsamples = 0
        for year, sample_map in year_sample_map:
            merged_dir = "/nfs-7/userdata/{}/tupler_babies/merged/FT/{}/output/year_{}".format(os.getenv("USER"),tag,year)
            for dsname,shortname in sample_map:
                if "/hadoop/" in dsname:
                    sample = DirectorySample(
                        dataset=dsname.split("|",1)[0].strip(),
                        location=dsname.split("|",1)[1].strip(),
                        )
                else:
                    sample = SNTSample(
                        dataset=dsname,
                        exclude_tag_pattern="CMS4_V08-*", # ignore new samples by sicheng for 2016
                        tag=tag_match,
                        # tag="CMS4_V09-04-13", # if not specified, get latest tag
                        )
                skip_tail = "/SMS" in dsname
                # skip_tail = False
                task = CondorTask(
                        sample = sample,
                        files_per_output = split_func(dsname),
                        output_name = "output.root",
                        tag = tag,
                        min_completion_fraction = 0.90 if skip_tail else 1.0,
                        condor_submit_params = {
                            # "sites":"T2_US_UCSD,UCSB",  # I/O is hella faster
                            "sites":"T2_US_UCSD",  # I/O is hella faster
                            # "sites":"UAF",  # I/O is hella faster
                            "classads": [ 
                                ["metis_extraargs",extra_args],
                                ["JobBatchName","FT_{}_{}".format(year,shortname)],
                                ],
                            # "classads": [ ["metis_extraargs","--ignorebadfiles"], ],
                            },
                        cmssw_version = "CMSSW_9_4_9",
                        scram_arch = "slc6_amd64_gcc630",
                        input_executable = "inputs/condor_executable_metis.sh",
                        tarfile = "inputs/package.tar.xz",
                        special_dir = "FTbabies/{}/".format(tag),
                        recopy_inputs = False,
                )
                merge_task = LocalMergeTask(
                        input_filenames=task.get_outputs(),
                        output_filename="{}/{}.root".format(merged_dir,shortname),
                        ignore_bad = skip_tail,
                        )
                nsamples += 1
                if not task.complete():
                    task.additional_input_files = ["/home/users/namin/2018/fourtop/all/FTAnalysis/babymaking/batch/condor_chirp"]
                    task.process()
                else:
                    if not merge_task.complete():
                        merge_task.process()
                    else:
                        ndone += 1
                        
                total_summary[dsname.split("|")[0]] = task.get_task_summary()

        StatsParser(data=total_summary, webdir="~/public_html/dump/FTmetis/").do(
                custom_event_rate_parser=get_event_rate,
                )

        if ndone == nsamples:
            print "All done!"
            sys.exit()

        time.sleep(600 if i < 10 else 1.5*60*60)

