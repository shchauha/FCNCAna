#!/usr/bin/env python
import os
import sys
pwd = os.getcwd()
sys.path.insert(1,pwd+'/..//utils')
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
# procs_to_consider = []
procs_to_consider = [ # FIXME
    #"fakes",
    #"data"
    ]


basedirs ={

    #2016: "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31/output/year_2016_94x/",
    #2017: "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31/output/year_2017/",
    #2018: "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.31//output/year_2018/",

    #2016: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.24/year_2016/",

    2016: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2016_94x/",
    2017: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2017/",
    2018: "/nfs-7/userdata/shchauha-2/tupler_babies/merged/FT/v3.31/year_2018/",

    }

outputdir = "outputs_v3p31_Sep7" #

options = {
    
    # for SS
    2016: "  doSS Data2016 new2016FRBins quiet ",
    2017: "  doSS Data2017 quiet ",
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
        "dy": make_objs([
                basedirs[2016]+"DY_high.root",
                basedirs[2016]+"DY_low.root",
                ], options=options[2016]),
        "wjets": make_objs(basedirs[2016]+"WJets_HT*.root", options=options[2016]),
        "fakes": make_objs(basedirs[2016]+"Data*.root", options=options[2016]+" doFakes"),
        "flips": make_objs(basedirs[2016]+"Data*.root", options=options[2016]+" doFlips"),
        "ttz": make_objs([basedirs[2016]+"TTZnlo.root",
                          basedirs[2016]+"TTZLOW.root",]
                         , options=options[2016]),
        
        "wz": make_objs([
                basedirs[2016]+"WZ.root",
                ],options=options[2016]),
        #        "singletop": make_objs([
        #                basedirs[2016]+"ST1.root",
        #                basedirs[2016]+"ST2.root",
        #                ],options=options[2016]),
        "ttdl": make_objs([
                basedirs[2016]+"TTDL.root",             
                ], options=options[2016]),
        
        "ttsl": make_objs([        
                basedirs[2016]+"TTSLtop.root",
                basedirs[2016]+"TTSLtopbar.root",
                ], options=options[2016]),                                
        
        "fcnc": make_objs(basedirs[2016]+"FCNC*.root", options=options[2016]),
        },
    2017: {
        "data": make_objs(basedirs[2017]+"Data*.root", options=options[2017]),
        "ttw": make_objs(basedirs[2017]+"TTWnlo.root", options=options[2017]),
        "tth": make_objs(basedirs[2017]+"TTHtoNonBB.root", options=options[2017]),
        "dy": make_objs([
                basedirs[2017]+"DY_high.root",
                basedirs[2017]+"DY_low.root",
                ], options=options[2017]),
        "wjets": make_objs(basedirs[2017]+"WJets_HT*.root", options=options[2017]),
        
        "fakes": make_objs(basedirs[2017]+"Data*.root", options=options[2017]+" doFakes"),
        "flips": make_objs(basedirs[2017]+"Data*.root", options=options[2017]+" doFlips"),
        "ttz": make_objs([basedirs[2017]+"TTZnlo.root",
                          basedirs[2017]+"TTZLOW.root",]
                         , options=options[2017]),
        
        "wz": make_objs([
                basedirs[2017]+"WZ.root",
                ],options=options[2017]),
        #"singletop": make_objs([
        #        basedirs[2017]+"ST1.root",
        #        basedirs[2017]+"ST2.root",
        #        ],options=options[2017]),
        "ttdl": make_objs([
                basedirs[2017]+"TTDL.root",             
                ], options=options[2017]),
        
        "ttsl": make_objs([        
                basedirs[2017]+"TTSLtop.root",
                basedirs[2017]+"TTSLtopbar.root",
                ], options=options[2017]),                        


        "fcnc": make_objs(basedirs[2017]+"FCNC*.root", options=options[2017]),

        },
    2018: {
        "data": make_objs([basedirs[2018]+"ReRecoData*.root",
                           basedirs[2018]+"Data*Dv2.root",
                           ], options=options[2018]),
        "ttw": make_objs(basedirs[2018]+"TTWnlo.root", options=options[2018]),
        "tth": make_objs(basedirs[2018]+"TTHtoNonBB.root", options=options[2018]),
        "dy": make_objs([
                basedirs[2018]+"DY_high.root",
                basedirs[2018]+"DY_low.root",
                ], options=options[2018]),
        "wjets": make_objs(basedirs[2018]+"WJets*.root", options=options[2018]),
        "fakes": make_objs([basedirs[2018]+"ReRecoData*.root",
                           basedirs[2018]+"Data*Dv2.root",
                           ], options=options[2018]+" doFakes"),        

        "flips": make_objs([basedirs[2018]+"ReRecoData*.root",
                           basedirs[2018]+"Data*Dv2.root",
                           ], options=options[2018]+" doFlips"),        
        "ttz": make_objs([basedirs[2018]+"TTZnlo.root",
                          basedirs[2018]+"TTZLOW.root",]
                         , options=options[2018]),
        
        "wz": make_objs([
                basedirs[2018]+"WZ.root",
                ],options=options[2018]),
        #"singletop": make_objs([
        #        basedirs[2018]+"ST1.root",
        #        basedirs[2018]+"ST2.root",
        #        ],options=options[2018]),
        "ttdl": make_objs([
                basedirs[2018]+"TTDL.root",             
                ], options=options[2018]),
        
        "ttsl": make_objs([        
                basedirs[2018]+"TTSLtop.root",
                basedirs[2018]+"TTSLtopbar.root",
                ], options=options[2018]),                        

        "fcnc": make_objs(basedirs[2017]+"FCNC*.root", options=options[2018]),
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
