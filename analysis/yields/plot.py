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

path = "outputs2016"
signalname = "fcnc"

files = []
for r, d, f in os.walk(path):
    for file in f:
        if '.root' in file:
            files.append(os.path.join(r, file))

log_scale =[False, True]
log_string =""
regions = ["sshh","ssbr"]


for region in regions:
    plotname = [
        [region+"_sr_in","SR"],
        [region+"_njets_in","No of jets"]    
        ]    
    for plot in range(len(plotname)):
        print plotname[plot][0]            
        bgs = []
        sigs = []
        data = []
        for f in files:
            print f
            if "data" in f:
                data.append(Hist1D(uproot.open(f)[plotname[plot][0]]))
            elif signalname in f:
                sigs.append(Hist1D(uproot.open(f)[plotname[plot][0]],label=signalname))
            else:
                bgs.append(Hist1D(uproot.open(f)[plotname[plot][0]],label=f.replace(".root","").replace(path+"/histos_","")))
        bgs = sorted(bgs, key=lambda bg: bg.get_integral())       
        #print bgs, data, sigs
        for log in log_scale:
            if log is True:
                log_string = "_log"
            else:
                log_string = ""

            plot_stack(bgs=bgs, data = "",title="", xlabel=plotname[plot][1], ylabel="Events", filename=path+"/"+plotname[plot][0]+log_string+".png",
                       cms_type = "Preliminary",
                       lumi = 35.9,
                       ratio = sum(sigs).divide(sum(bgs)),
                       ratio_range=[0,1.5],
                       sigs=sigs,
                       do_log=log,    
                       mpl_ratio_params={"label":"Sig/Pred."},
                       
                       )
        
        
