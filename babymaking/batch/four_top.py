import os
import sys

d_ds2name = {
        "backgrounds": {

            "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM":"TTBAR_PH",

            # XXX note, add scale1fbs when adding new samples, and also copy the package file

            # "/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTSL",
            # "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTBAR_PH",
            # "/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "W3Jets",
            # "/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "W4Jets",

            # "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM":"DY_high",
            # "/TTTT_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM":"TTTTnew",
            # "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM":"TTBAR_PH",

            # "/TTTT_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_pilot_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTTTnew",
            # "/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTZnlo",
            # "/TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTZLOW",
            # "/ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM": "TTHtoNonBB",
            # "/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTWnlo",
            # "/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM": "DY_low",
            # "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "DY_high",
            # "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM": "DY_highext1",
            # "/WW_TuneCP5_13TeV-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM": "WW",
            # "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "ZZ",
            # "/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "WZZ",
            # "/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "ZZZ",
            # "/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "WWZ",
            # "/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "WWW",
            # "/TTTJ_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTTJ",
            # "/TTTW_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTTW",
            # "/TTHH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTHH",
            # "/TTWW_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM": "TTWW",
            # "/TTZZ_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTZZ",
            # "/TTWZ_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTWZ",
            # "/TTWH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTWH",
            # "/TTZH_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TTZH",
            # "/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "TZQ",
            # "/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "WZG",
            # "/WW_DoubleScattering_13TeV-pythia8_TuneCP5/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "WWDPS",
            # "/ST_tW_top_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "STtop",
            # "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM": "STantitop",
            # "/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v2/MINIAODSIM": "WZ",
            # # "/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM": "GGHtoZZto4L",
            # # "/ST_tWll_5f_LO_13TeV-MadGraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM": "TWZ",
            # # "/ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM": "ZG",
            # # "/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM": "TTG",
            # # "/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM": "QQWW",
            # # "/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM": "TGext",
            # # "/WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM": "WGToLNuGext",
            # # "/VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM": "VHtoNonBB",
            # # "/WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM": "WWG",
            # # "/TTGamma_SingleLeptFromT_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM": "TTGsinglelep",
            # # "/TTGamma_SingleLeptFromTbar_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM": "TTGsinglelepbar",
            # # "/TTGamma_Dilept_TuneCUETP8M2T4_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM": "TTGdilep",

            # # note: 92X
            # "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v1/MINIAODSIM": "WJets",
            # "/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM": "W4Jets",
            # "/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v3/MINIAODSIM": "W3Jets",
            # "/TTToSemiLeptonic_mtop175p5_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM": "TTSLm175",
            # "/TTToSemiLeptonic_mtop171p5_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM": "TTSLm171",
            # "/TTToSemiLeptonic_mtop173p5_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM": "TTSLm173",



            },
        "data": {

                # "/MuonEG/Run2016C-17Jul2018-v1/MINIAOD": "DataMuonEGC",

                # "/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD": "DataDoubleEGB",
                # "/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD": "DataDoubleEGC",
                # "/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD": "DataDoubleEGD",
                # "/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD": "DataDoubleEGE",
                # "/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD": "DataDoubleEGF",
                # "/MuonEG/Run2017B-31Mar2018-v1/MINIAOD": "DataMuonEGB",
                # "/MuonEG/Run2017C-31Mar2018-v1/MINIAOD": "DataMuonEGC",
                # "/MuonEG/Run2017D-31Mar2018-v1/MINIAOD": "DataMuonEGD",
                # "/MuonEG/Run2017E-31Mar2018-v1/MINIAOD": "DataMuonEGE",
                # "/MuonEG/Run2017F-31Mar2018-v1/MINIAOD": "DataMuonEGF",
                # "/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD": "DataDoubleMuonB",
                # "/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD": "DataDoubleMuonC",
                # "/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD": "DataDoubleMuonD",
                # "/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD": "DataDoubleMuonE",
                # "/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD": "DataDoubleMuonF",

            },

        }

def dataset_to_shortname(ds):
    for typ_values in d_ds2name.values():
        for dname in typ_values:
            if dname.split("|")[0] == ds:
                return typ_values[dname]
    return "unknown"

inputs_path = os.path.abspath(__file__).rsplit("/",1)[0]+"/inputs/"
sweepRoot_scripts = [inputs_path+"sweeproot.sh", inputs_path+"sweeproot_macro.C"]
merging_scripts = [inputs_path+"merge_script.sh", inputs_path+"merge_macro.C"]
baby_merged_dir = "/nfs-7/userdata/${USER}/tupler_babies/merged/"

# tag = "v1.00_94x_baseline_v1"
# tag = "v1.00_94x_baseline_dyonly_v1"
# tag = "v1.00_94x_baseline_withos_v3"
# tag = "v1.00_94x_baseline_withos_v3_2017egid"
# tag = "v1.02_94x_withos_relaxhits"
# tag = "v1.02_94x_withos_relaxhitsmore"
# tag = "v1.03_os_nmiss1_dileptrig_v1"
# tag = "v1.04_v1" # new 2017 ELE MVA, new samples, new effective areas
# tag = "v1.04_v2" # fix 2017 ELE MVA, more samples
# tag = "v1.04_v2_dyeedata" # same as above but only DY and DoubleEG data
# tag = "v1.05_v1" # updated multiIso WPs and rhos
# tag = "v1.06_v1" # corrected conecorrpt calculation
# tag = "v1.06_v2" # same as v1 but with modified met (https://indico.cern.ch/event/718554/contributions/3028859/attachments/1668190/2675164/METin2017_June14_2018_Vesterbacka.pdf)
# tag = "v1.06_v2_oldiso" # only different wrt v1.06_v2 is I've temporarily(!!!) reverted to the 2016 multi isolation cuts -- everything else is identical -- only ttbar,Wjets
# tag = "v1.06_v2_oldisofix" # diff wrt previous line is I corrected conecorr pt constants to match the old WPs as well (very small -- % level -- difference I hope)
# tag = "v1.06_v3" # same as v2 but testing tttt,dy,tt samples with fixed nmiss value
# tag = "v1.06_v2_sleptrigv2" # same as v2 but added branch for iso single lep trigs
# tag = "v1.06_v2_sleptrigv2" # previous line ran over new CMS4 which missed some files, so try 2
# tag = "v2.00_94x_muegdata" # trying out run2 babymaker with 94x rereco of 2016C MuEG
tag = "v2.00_2017_insitu_v2" # 2017 MC ttbar with insitu hyp_classes, added needed tight MVA branch
package = inputs_path+"package.tar.gz"
executable = inputs_path+"condor_executable.sh"
dashboard_name = "AutoTwopler_FTbabies"

if not os.path.isfile(package):
    raise Exception("Woah! Did you run make_tar.sh? I don't see the package.tar.gz file in %s" % inputs_path)
elif ((os.path.getmtime(package) < os.path.getmtime("main.exe")) or (os.path.getmtime(package) < os.path.getmtime("helper_babymaker.cc"))):
    raise Exception("Woah! Did you run make_job_inputs.sh? It looks like main.exe or helper_babymaker.cc might be newer than the tarball in %s" % inputs_path)



