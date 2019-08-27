from tqdm import tqdm
import json
import os
import numpy as np
import itertools
import uproot
import re
#from analysis.utils.plotting_utils import write_table
import sys
sys.path.insert(0,'/home/users/{}/.local/lib/python2.7/site-packages/'.format(os.getenv("USER")))
from matplottery.plotter import plot_stack
from matplottery.utils import Hist1D, MET_LATEX, binomial_obs_z

#dirname = "cutandcountmc"
#dirname = "outputs2016"
#dirname = "v3.24_test"
#dirname = "v4.2"
dirname = "v4.2_data"


signalname = "fcnc"
files = []

for r, d, f in os.walk(dirname):
    for file in f:
        if file.endswith(".root"):
            files.append(file)
print files
log_scale =[False, True]
log_string =""
regions = ["ssbr","ss0b2j","ss1b2j","ss2b2j"]
#regions = ["ss1b2j_s_met50"]

for region in regions:
    print "working on ",region
    plotname = [
        #[region+"_sr_in","SR"],
        [region+"_njets_in","No of jets"],
        [region+"_nbtags_in","No of b jets"],
        [region+"_met_in","MET"],
        [region+"_ht_in","HT"],
        [region+"_mtmin_in","MTmin"],
        [region+"_mll_in","M_{ll} [GeV]"],
        [region+"_dphil1met_in","dphil1met"],
        [region+"_dphil2met_in","dphil2met"],
        [region+"_dphimetj1_in","dphimetj1"],
        [region+"_ptrel1_in","ptrel1"],
        [region+"_ptrel2_in","ptrel2"],
        [region+"_ptratio1_in","ptratio1"],
        [region+"_ptratio2_in","ptratio2"],
        [region+"_miniiso1_in","miniiso1"],
        [region+"_miniiso2_in","miniiso2"],
        [region+"_dphil1l2_in","dphil1l2"],
        [region+"_bdt_in","bdt"],
        [region+"_ptj1_in","1st jet pt"],
        [region+"_ptj2_in","2nd jet pt"],
        [region+"_ptj3_in","3rd jet pt"],
        [region+"_ptj4_in","4th jet pt"],
        [region+"_ptbt1_in","1st bjet pt"],
        [region+"_ptbt2_in","2nd bjet pt"],
        [region+"_ptbt3_in","3rd bjet pt"],
        [region+"_ptbt4_in","4th bjet pt"],

        ]    
    for plot in range(len(plotname)):
        #print plotname[plot][0]            
        bgs = []
        sigs = []
        data = []
        for f in files:
            #print f
            if "data" in f:
                data.append(Hist1D(uproot.open(dirname+"/"+f)[plotname[plot][0]]))
            elif signalname in f:
                sigs.append(Hist1D(uproot.open(dirname+"/"+f)[plotname[plot][0]],label=signalname))
            else:
                bgs.append(Hist1D(uproot.open(dirname+"/"+f)[plotname[plot][0]],label=f.replace(".root","").replace("histos_","")))
        bgs = sorted(bgs, key=lambda bg: bg.get_integral())       
        #print bgs, data, sigs
        for log in log_scale:
            if log is True:
                log_string = "_log"
            else:
                log_string = ""

            plot_stack(bgs=bgs,
                       #data = data[0],
                       title="",
                       xlabel=plotname[plot][1], 
                       ylabel="Events", 
                       filename=dirname+"/"+plotname[plot][0]+log_string+".png",
                       cms_type = "Preliminary",
                       #lumi = 140,
                       lumi = 35.9,
                       #ratio = data[0].divide(sum(bgs)),
                       ratio = sum(sigs).divide(sum(bgs)),
                       ratio_range=[0,2],
                       sigs=sigs,
                       do_log=log,    
                       mpl_ratio_params={"label":"Sig/Bkgd"},                       
                       )
        
        del bgs
        del sigs
        del data
