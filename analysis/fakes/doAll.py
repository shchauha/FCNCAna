#!/usr/bin/env python
import os
import sys
sys.path.insert(1,os.getcwd()+'/..//utils')
import pyrun as pyrun
import argparse
import fnmatch
import operator
import glob
import ast
import time
import ROOT as r
r.gROOT.SetBatch()
r.gROOT.ProcessLine(".L ../misc/class_files/v8.02/SS.cc+")
r.gROOT.ProcessLine(".L ../../common/CORE/Tools/dorky/dorky.cc+")
r.gROOT.ProcessLine(".L ScanChain.C+")
years_to_consider = [ # FIXME
    2016,
    2017,
    2018,
]
#procs_to_consider = []
procs_to_consider = [ # FIXME    
    #"fcnc_hut",
    #'ST_hut'
    "data",
    "truthfakes_ttsl",
    "truthfakes_wjets",
    "fakes_ttsl",
    "fakes_wjets"
    
]

basedirs ={
    #2016: "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31/output/year_2016_94x/",
    #2017: "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31/output/year_2017/",
    #2018: "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31//output/year_2018/",
    #2016: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.24/year_2016/",
    
    #2016: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2016_94x/",
    #2017: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2017/",
    #2018: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2018/",
    2016: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2016_94x/skim/",
    2017: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2017/skim/",
    2018: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2018/skim/",
    
}

#outputdir = "outputs_v3p31_SR_FCNC_Jan8_MLtraining" #
#outputdir = "outputs_v3p31_ttz" #
#outputdir = "outputs_v3p31_fake_closure"
outputdir = "outputs_v3p31_apr27"
 
options = {    
    # for SS
    #2016: "  doSS Data2016 new2016FRBins ReadBDT quiet ",
    #2017: "  doSS Data2017 ReadBDT quiet ",
    #2018: "  doSS Data2018 ReadBDT quiet ",

    2016: "  doSS Data2016 new2016FRBins quiet ",
    2017: "  doSS Data2017  quiet ",
    2018: "  doSS Data2018 quiet ",
    }

def make_objs(fpatts=[],options="",treename="t"):
    if type(fpatts) == str: fpatts = [fpatts]
    ch = r.TChain(treename)
    for fpatt in fpatts:
        ch.Add(fpatt)
    return {"ch": ch, "options": options}

chs = {
    2016: {
        "data": make_objs(basedirs[2016]+"Data*.root", options=options[2016]),
        "ttw": make_objs(basedirs[2016]+"TTWnlo.root", options=options[2016]),
        "tth": make_objs(basedirs[2016]+"TTHtoNonBB.root", options=options[2016]),
        "ttz": make_objs([basedirs[2016]+"TTZnlo.root",
                          basedirs[2016]+"TTZLOW.root",]
                         , options=options[2016]),
        "fakes": make_objs([
                basedirs[2016]+"Data*.root",
                basedirs[2016]+"TTWnlo.root",
                basedirs[2016]+"TTZnlo.root",
                basedirs[2016]+"TTHtoNonBB.root",
                ] , options=options[2016]+" doFakes doData"),

        "fakes_exp": make_objs([
                basedirs[2016]+"TTSLtop.root",
                basedirs[2016]+"TTSLtopbar.root",
                basedirs[2016]+"WJets_HT*.root",
                ] , options=options[2016]+" doFakes doExp"),

        "fakes_mc": make_objs([
                basedirs[2016]+"TTSLtop.root",
                basedirs[2016]+"TTSLtopbar.root",
                basedirs[2016]+"WJets_HT*.root",
                ] , options=options[2016]+ " doTruthFake"),
        
        "truthfakes_ttsl": make_objs([
            basedirs[2016]+"TTSLtop.root",
            basedirs[2016]+"TTSLtopbar.root",
            #basedirs[2016]+"WJets_HT*.root",
        ] , options=options[2016]+"  doTruthFake"),

        "truthfakes_wjets": make_objs([
            #basedirs[2016]+"TTSLtop.root",
            #basedirs[2016]+"TTSLtopbar.root",
            basedirs[2016]+"WJets_HT*.root",
        ] , options=options[2016]+" doTruthFake"),
    
        "fakes_ttsl": make_objs([
            basedirs[2016]+"TTSLtop.root",
            basedirs[2016]+"TTSLtopbar.root",
            #basedirs[2016]+"WJets_HT*.root",
        ] , options=options[2016]+" doFakes doTruthFake"),

        "fakes_wjets": make_objs([
            #basedirs[2016]+"TTSLtop.root",
            #basedirs[2016]+"TTSLtopbar.root",
            basedirs[2016]+"WJets_HT*.root",
        ] , options=options[2016]+" doFakes doTruthFake"),
    

        "flips": make_objs(basedirs[2016]+"Data*.root", options=options[2016]+" doFlips"),        
        "flips_mc": make_objs([
                basedirs[2016]+"TTDL.root",
                basedirs[2016]+"DY_high.root",
                basedirs[2016]+"DY_low.root",
                ] ,options=options[2016]),
        "fakes_mc_ml": make_objs([
                basedirs[2016]+"TTDL.root",
                basedirs[2016]+"DY_high.root",
                basedirs[2016]+"DY_low.root",
                ] ,options=options[2016]+ " doTruthFake"),
        "xg": make_objs([
                basedirs[2016]+"TGext.root",
                basedirs[2016]+"TTGdilep.root",
                basedirs[2016]+"TTGsinglelepbar.root",
                basedirs[2016]+"TTGsinglelep.root",
                basedirs[2016]+"WGToLNuGext.root",
                basedirs[2016]+"ZG.root",
                ],options=options[2016] + " doXgamma "),
        "ttvv": make_objs([
                basedirs[2016]+"TTHH.root",
                basedirs[2016]+"TTWH.root",
                basedirs[2016]+"TTWW.root",
                basedirs[2016]+"TTWZ.root",
                basedirs[2016]+"TTZH.root",
                basedirs[2016]+"TTZZ.root",
                ],options=options[2016]),
        "rares": make_objs([
                basedirs[2016]+"GGHtoZZto4L.root",
                basedirs[2016]+"QQWW.root",
                basedirs[2016]+"TWZ.root",
                basedirs[2016]+"TZQ.root",
                basedirs[2016]+"VHtoNonBB.root",
                basedirs[2016]+"WWDPS.root",
                basedirs[2016]+"WWW.root",
                basedirs[2016]+"WWZ.root",
                basedirs[2016]+"WZ.root",
                basedirs[2016]+"WZG.root",
                basedirs[2016]+"WWG.root",
                basedirs[2016]+"WZZ.root",
                basedirs[2016]+"ZZ.root",
                basedirs[2016]+"ZZZ.root",
                basedirs[2016]+"TTTJ.root",
                basedirs[2016]+"TTTW.root",
                ],options=options[2016]),
        "fcnc_hut": make_objs(basedirs[2016]+"FCNC_hut*.root", options=options[2016]),
        "fcnc_hct": make_objs(basedirs[2016]+"FCNC_hct*.root", options=options[2016]),
        "ST_hut": make_objs("/nfs-7/userdata/shchauha/tupler_babies/merged/FT/v3.31/output/year_2016/ST_hut_top.root", options=options[2016]),
        },
    2017: {
        "data": make_objs(basedirs[2017]+"Data*.root", options=options[2017]),
        "ttw": make_objs(basedirs[2017]+"TTWnlo.root", options=options[2017]),
        "tth": make_objs(basedirs[2017]+"TTHtoNonBB.root", options=options[2017]),
        "ttz": make_objs([basedirs[2017]+"TTZnlo.root",
                          basedirs[2017]+"TTZLOW.root",]
                         , options=options[2017]),
        "fakes": make_objs([
                basedirs[2017]+"Data*.root",
                basedirs[2017]+"TTWnlo.root",
                basedirs[2017]+"TTZnlo.root",
                basedirs[2017]+"TTHtoNonBB.root",
                ] , options=options[2017]+" doFakes doData"),

        "fakes_exp": make_objs([
                basedirs[2017]+"TTSLtop.root",
                basedirs[2017]+"TTSLtopbar.root",
                basedirs[2017]+"WJets_HT*.root",
                
                ] , options=options[2017]+" doFakes doExp"),

        "fakes_mc": make_objs([
                basedirs[2017]+"TTSLtop.root",
                basedirs[2017]+"TTSLtopbar.root",
                basedirs[2017]+"WJets_HT*.root",
                ] , options=options[2017]+ " doTruthFake"),

        "truthfakes_ttsl": make_objs([
            basedirs[2017]+"TTSLtop.root",
            basedirs[2017]+"TTSLtopbar.root",
            #basedirs[2017]+"WJets_HT*.root",
        ] , options=options[2017]+"  doTruthFake"),

        "truthfakes_wjets": make_objs([
            #basedirs[2017]+"TTSLtop.root",
            #basedirs[2017]+"TTSLtopbar.root",
            basedirs[2017]+"WJets_HT*.root",
        ] , options=options[2017]+" doTruthFake"),
    
        "fakes_ttsl": make_objs([
            basedirs[2017]+"TTSLtop.root",
            basedirs[2017]+"TTSLtopbar.root",
            #basedirs[2017]+"WJets_HT*.root",
        ] , options=options[2017]+" doFakes doTruthFake"),

        "fakes_wjets": make_objs([
            #basedirs[2017]+"TTSLtop.root",
            #basedirs[2017]+"TTSLtopbar.root",
            basedirs[2017]+"WJets_HT*.root",
        ] , options=options[2017]+" doFakes doTruthFake"),
    

        "flips": make_objs(basedirs[2017]+"Data*.root", options=options[2017]+" doFlips"),        
        "flips_mc": make_objs([
                basedirs[2017]+"TTDL.root",
                basedirs[2017]+"DY_high.root",
                basedirs[2017]+"DY_low.root",
                ] ,options=options[2017]),
        "fakes_mc_ml": make_objs([
                basedirs[2017]+"TTDL.root",
                basedirs[2017]+"DY_high.root",
                basedirs[2017]+"DY_low.root",
                ] ,options=options[2017]+ " doTruthFake"),
        "xg": make_objs([
                basedirs[2017]+"TGext.root",
                basedirs[2017]+"TTGdilep.root",
                basedirs[2017]+"TTGsinglelepbar.root",
                basedirs[2017]+"TTGsinglelep.root",
                basedirs[2017]+"WGToLNuGext.root",
                basedirs[2017]+"ZG.root",
                ],options=options[2017] + " doXgamma "),
        "ttvv": make_objs([
                basedirs[2017]+"TTHH.root",
                basedirs[2017]+"TTWH.root",
                basedirs[2017]+"TTWW.root",
                basedirs[2017]+"TTWZ.root",
                basedirs[2017]+"TTZH.root",
                basedirs[2017]+"TTZZ.root",
                ],options=options[2017]),
        "rares": make_objs([
                basedirs[2017]+"GGHtoZZto4L.root",
                basedirs[2017]+"QQWW.root",
                basedirs[2017]+"TWZ.root",
                basedirs[2017]+"TZQ.root",
                basedirs[2017]+"VHtoNonBB.root",
                basedirs[2017]+"WWDPS.root",
                basedirs[2017]+"WWW.root",
                basedirs[2017]+"WWZ.root",
                basedirs[2017]+"WZ.root",
                basedirs[2017]+"WZG.root",
                basedirs[2017]+"WWG.root",
                basedirs[2017]+"WZZ.root",
                basedirs[2017]+"ZZ.root",
                basedirs[2017]+"ZZZ.root",
                basedirs[2017]+"TTTJ.root",
                basedirs[2017]+"TTTW.root",
                ],options=options[2017]),
        #"fcnc": make_objs(basedirs[2017]+"FCNC*.root", options=options[2017]),
        "fcnc_hut": make_objs(basedirs[2017]+"FCNC_hut*tauDecay.root", options=options[2017]), # samples for tauDecay
        "fcnc_hct": make_objs(basedirs[2017]+"FCNC_hct*tauDecay.root", options=options[2017]), # samples for tauDecay
        "ST_hut": make_objs("/nfs-7/userdata/shchauha/tupler_babies/merged/FT/v3.31/output/year_2017/ST_hut_top.root", options=options[2017]),
        },
    2018: {
        "data": make_objs([basedirs[2018]+"ReRecoData*.root",
                           basedirs[2018]+"Data*Dv2.root",
                           ], options=options[2018]),
        "ttw": make_objs(basedirs[2018]+"TTWnlo.root", options=options[2018]),
        "tth": make_objs(basedirs[2018]+"TTHtoNonBB.root", options=options[2018]),
        "ttz": make_objs([basedirs[2018]+"TTZnlo.root",
                          basedirs[2018]+"TTZLOW.root",]
                         , options=options[2018]),
        "fakes": make_objs([
                basedirs[2018]+"ReRecoData*.root",
                basedirs[2018]+"Data*Dv2.root",
                basedirs[2018]+"TTWnlo.root",
                basedirs[2018]+"TTZnlo.root",
                basedirs[2018]+"TTHtoNonBB.root",
                ] , options=options[2018]+" doFakes doData"),

        "fakes_exp": make_objs([
                basedirs[2018]+"TTSLtop.root",
                basedirs[2018]+"TTSLtopbar.root",
                basedirs[2018]+"WJets_HT*.root",

                ] , options=options[2018]+" doFakes doExp"),

        "fakes_mc": make_objs([
                basedirs[2018]+"TTSLtop.root",
                basedirs[2018]+"TTSLtopbar.root",
                basedirs[2018]+"WJets_HT*.root",
                ] , options=options[2018]+ " doTruthFake"),

        "truthfakes_ttsl": make_objs([
            basedirs[2018]+"TTSLtop.root",
            basedirs[2018]+"TTSLtopbar.root",
            #basedirs[2017]+"WJets_HT*.root",
        ] , options=options[2018]+"  doTruthFake"),

        "truthfakes_wjets": make_objs([
            #basedirs[2017]+"TTSLtop.root",
            #basedirs[2017]+"TTSLtopbar.root",
            basedirs[2018]+"WJets_HT*.root",
        ] , options=options[2018]+" doTruthFake"),
    
        "fakes_ttsl": make_objs([
            basedirs[2018]+"TTSLtop.root",
            basedirs[2018]+"TTSLtopbar.root",
            #basedirs[2017]+"WJets_HT*.root",
        ] , options=options[2018]+" doFakes doTruthFake"),

        "fakes_wjets": make_objs([
            #basedirs[2017]+"TTSLtop.root",
            #basedirs[2017]+"TTSLtopbar.root",
            basedirs[2018]+"WJets_HT*.root",
        ] , options=options[2018]+" doFakes doTruthFake"),
    


        "flips": make_objs([                
                basedirs[2018]+"ReRecoData*.root",                
                basedirs[2018]+"Data*Dv2.root",   
                ], options=options[2018]+" doFlips"),        
        "flips_mc": make_objs([
                basedirs[2018]+"TTDL.root",
                basedirs[2018]+"DY_high.root",
                basedirs[2018]+"DY_low.root",
                ] ,options=options[2018]),
        "fakes_mc_ml": make_objs([
                basedirs[2018]+"TTDL.root",
                basedirs[2018]+"DY_high.root",
                basedirs[2018]+"DY_low.root",
                ] ,options=options[2018]+ " doTruthFake"),        
#        "wz": make_objs([
#                basedirs[2018]+"WZ.root",
#                ],options=options[2018]),        
        "xg": make_objs([
                basedirs[2018]+"TGext.root",
                basedirs[2018]+"TTGdilep.root",
                basedirs[2018]+"TTGsinglelepbar.root",
                basedirs[2018]+"TTGsinglelep.root",
                basedirs[2018]+"WGToLNuGext.root",
                basedirs[2018]+"ZG.root",
                ],options=options[2018] + " doXgamma "),
        "ttvv": make_objs([
                basedirs[2018]+"TTHH.root",
                basedirs[2018]+"TTWH.root",
                basedirs[2018]+"TTWW.root",
                basedirs[2018]+"TTWZ.root",
                basedirs[2018]+"TTZH.root",
                basedirs[2018]+"TTZZ.root",
                ],options=options[2018]),
        "rares": make_objs([
                basedirs[2018]+"GGHtoZZto4L.root",
                basedirs[2018]+"QQWW.root",
                basedirs[2018]+"TWZ.root",
                basedirs[2018]+"TZQ.root",
                basedirs[2018]+"VHtoNonBB.root",
                basedirs[2018]+"WWDPS.root",
                basedirs[2018]+"WWW.root",
                basedirs[2018]+"WWZ.root",
                basedirs[2018]+"WZ.root",
                basedirs[2018]+"WZG.root",
                basedirs[2018]+"WWG.root",
                basedirs[2018]+"WZZ.root",
                basedirs[2018]+"ZZ.root",
                basedirs[2018]+"ZZZ.root",
                basedirs[2018]+"TTTJ.root",
                basedirs[2018]+"TTTW.root",
                ],options=options[2018]),        
        "fcnc_hut": make_objs(basedirs[2018]+"FCNC_hut*tauDecay.root", options=options[2018]), # samples for tauDecay
        "fcnc_hct": make_objs(basedirs[2018]+"FCNC_hct*tauDecay.root", options=options[2018]), # samples for tauDecay
        "ST_hut": make_objs("/nfs-7/userdata/shchauha/tupler_babies/merged/FT/v3.31/output/year_2017/ST_hut_top.root", options=options[2018]),
        }
    }

# Change chain titles to proc_year so that we output the right root file name
for year in chs:
    [obj["ch"].SetTitle("{}_{}".format(proc,year)) for proc,obj in chs[year].items()]

def run_chain((index,info)):
    ch, options, outputdir = info
    ret = r.ScanChain(ch, options, outputdir)
    return index, ret

to_run = []
for year in years_to_consider:
    for proc,obj in chs[year].items():
        if procs_to_consider and (proc not in procs_to_consider): continue
        to_run.append([obj["ch"], obj["options"], outputdir])

os.system("mkdir -p {}".format(outputdir))

os.nice(10)
runner = pyrun.Runner(nproc=min(len(to_run),25), func=run_chain, dot_type=2)
runner.add_args(to_run)
runner.run()
