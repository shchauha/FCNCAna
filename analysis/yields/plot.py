#!/usr/bin/env python
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
#dirname = "v3.24_training"
#dirname = "v4.2/mc"
#dirname = "v4.2_nonskim"
#dirname = "v4.2_data"
#dirname = "v4.2_skimfix2"
#dirname = "v3.31_2016_94x_140ifb"
#dirname = "v3.31_2016_94x_140ifb/data"
#dirname = "v3.31_2016_94x/mc"
#dirname = "v3.31_2016_94x/data"
#dirname = "v3.31_2016_94x_140ifb_plots/mc"
#dirname = "v3.31_2017/mc"
#dirname = "v3.31_2017/data"

#dirname = "v3.31_2016/osmc"
dirname = "outputs_v3p31_Sep7"
plots_dir = dirname+"/mc"
os.system("mkdir -p {}".format(plots_dir))
#proc_to_consider = ["ttw","tth","ttz","ttdl","ttsl","wz","singletop","wjets","dy"]
proc_to_consider = ["ttw","tth","ttz","ttdl","ttsl","wz","wjets","dy"]
#proc_to_consider = ["ttw","tth","ttz","fakes","flips","wz"]

signalname = "fcnc"
#signalname = "ttw"

year = "2016"
other_years = ["2017","2018"]
#other_years = []

files = {}
for y in [year]+other_years:
    files[y] = { }
    for proc in proc_to_consider+["data"]+[signalname]:
        # try:     
        ystr = str(y)
        # if any(x in proc for x in ["fs_","rpv_"]): ystr = "2016"          
        files[y][proc] = uproot.open(dirname+"/histos_"+proc+"_"+ystr+".root")
#dy= files["2016"]["dy"]    

#print files

log_scale =[False, True]
log_string =""
regions = [
    "sshh",
    "ssbr",
    "ss0b2j",
    "ss1b2j",
    "ss2b2j",
    "mlbr",
    "ml1b1j",
    "ml2b2j",
    "mlbronz",
    "ml1b1jonz",
    "ml2b2jonz",
    "mlbroffz",
    "ml1b1joffz",
    "ml2b2joffz",
    "osbr",
    "tl",
    "br",
    "susytl",
    "hhtl",
    #"lowmetonzor0b",
    #"ssbr2",
    #"ss1b2j2",
    #"ss2b2j2",
    #           "ss1b2jbtagM",
    #           "ss2b2jbtagM",
    #           "ss1b2jbtag25",
    #           "ss2b2jbtag25",
    #           "ss1b2jbtag25M",
    #           "ss2b2jbtag25M",
    #           "ss1b2jjet40",
    #           "ss2b2jjet40",
    #           "ss1b2jjet40btagM",
    #           "ss2b2jjet40btagM",
    #           "ss1b2jjet40btag25",
    #           "ss2b2jjet40btag25",
    #           "ss1b2jjet40btag25M",
    #           "ss2b2jjet40btag25M",
    #
    #
    ]

d_label_colors = {
    "dy" : [0.4, 0.6, 1.0],
    "fakes" :  [0.85, 0.85, 0.85],
    "ttsl" :  [0.85, 0.85, 0.85],
    "singletop": [1.0, 0.4, 0.0],
    # "fakes_mc" : (r"MC fakes", [0.85, 0.85, 0.85]),
    "flips" :  [0.4, 0.4, 0.4],
    "ttdl" :  [0.4, 0.4, 0.4],
    "rares" :  [1.0, 0.4, 1.0],
    "tth"   :  [0.4, 0.4, 0.6],
    "wz"  :  [1.0,0.8,0.0],
    "wjets"  : [1.0,0.6,0.0],
    "ttw"   :  [0.0, 0.4, 0.0],
    "ttz"   : [0.4, 0.8, 0.4],
    "xg"    :  [0.4, 0.0, 0.8],
    }

d_flat_systematics = {
    "fakes": 0.40,
    "flips": 0.2,
    "ttsl": 0.40,
    "wjets": 0.40,
    "singletop": 0.40,
    "ttdl": 0.2,
    "rares": 0.5,
    "ww": 0.3,
    "ttw": 0.3,
    "ttz": 0.3,
    "wz": 0.3,
    "tth": 0.3,
    "xg": 0.5,
    }


for region in regions:
    print "working on ",region
    plotname = [
        #[region+"_sr_in","SR"],
        [region+"_njets_in","No of jets"],
        [region+"_nbtags_in","No of b jets"],
        #        [region+"_njets40_in","No of jets pt 40"],
        #        [region+"_nbtags25_in","No of bjets pt 25"],
        #        [region+"_nbtags25M_in","No of bjets pt 25 medium"],
        #        [region+"_nbtagsM_in","No of bjets medium"],
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
        #        [region+"_bdt_in","bdt"],
        #
        [region+"_drl1l2_in","drl1l2"],
        [region+"_mindrl1j_in","mindrl1j"],
        [region+"_mindrl2j_in","mindrl2j"],
        [region+"_mindrl1bt_in","mindrl1bt"],
        [region+"_mindrl2bt_in","mindrl2bt"],
        [region+"_mt1_in","mt1"],        
        [region+"_mt2_in","mt2"],
        [region+"_l1dxy_in","l1dxy"],
        [region+"_l2dxy_in","l2dxy"],
        [region+"_l1dz_in","l1dz"],
        [region+"_l2dz_in","l2dz"],
        [region+"_mossf_in","mossf"],
        #        #[region+"__in",""],        
        [region+"_ptj1_in","1st jet pt"],
        [region+"_ptj2_in","2nd jet pt"],
        [region+"_ptj3_in","3rd jet pt"],
        #[region+"_ptj4_in","4th jet pt"],
        [region+"_ptbt1_in","1st bjet pt"],
        [region+"_ptbt2_in","2nd bjet pt"],
        #        [region+"_ptbt3_in","3rd bjet pt"],
        #        [region+"_ptbt4_in","4th bjet pt"],                
        [region+"_fwd_jetpt_in","fwd jet pt"],        
        #

        ]    
    for plot in range(len(plotname)):
        data = sum([Hist1D(files[y]["data"][plotname[plot][0]], label = "data") for y in files.keys()])
        sigs = [sum([Hist1D(files[y][signalname][plotname[plot][0]], label = signalname, color = d_label_colors.get(signalname) ) for y in files.keys()])]
        bgs=[sum([Hist1D(files[y][proc][plotname[plot][0]], label = proc, color = d_label_colors.get(proc)) for y in files.keys()])
             for proc in proc_to_consider
             ]        
        bgs = sorted(bgs, key=lambda bg: bg.get_integral())       

        #print "size of the bgs ", len(bgs)
        for bg in bgs:
            # add flat systematic to stat unc in quadrature
            bg._errors = np.hypot(bg._counts*d_flat_systematics.get(bg.get_attr("label"),0.),bg._errors)
        #print  sigs
        for log in log_scale:
            if log is True:
                log_string = "_log"
            else:
                log_string = ""
                
            lumi = 35.9     
            if year is "2017":
                lumi = 41.5            
            if year is "2018":
                lumi = 59.7
            outname = plots_dir+"/year"+year+"_"+plotname[plot][0]+log_string+".png"
            if len(other_years) > 0:
                outname = plots_dir+"/run2"+"_"+plotname[plot][0]+log_string+".png"            
                lumi = 137
            print year, outname, lumi, len(other_years)    

            plot_stack(bgs=bgs,
                       data = data,
                       title="",
                       xlabel=plotname[plot][1], 
                       ylabel="Events",                        
                       filename=outname,
                       cms_type = "Preliminary",                       
                       lumi = lumi,                       
                       ratio = data.divide(sum(bgs)),
                       #ratio = sum(sigs).divide(sum(bgs)),
                       ratio_range=[0,2],
                       sigs=sigs,
                       do_log=log,    
                       #mpl_ratio_params={"label":"Sig/Bkgd"},                       
                       mpl_ratio_params={"label":"Obs./Pred."},                       
                       do_bkg_syst=True,
                       )
        
        del bgs
        del sigs
        del data
