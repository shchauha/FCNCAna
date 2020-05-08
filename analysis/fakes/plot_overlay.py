from tqdm import tqdm
import json
import os
import sys
import numpy as np
import itertools
import uproot

sys.path.insert(1,os.getcwd()+'/../utils')
from plotting_utils import write_table
sys.path.insert(0,os.getenv("HOME")+'/.local/lib/python2.7/site-packages/')
from matplottery.plotter import plot_stack
from matplottery.utils import Hist1D, MET_LATEX
np.set_printoptions(linewidth=200)

plotdata = False
#plotdata = True
useddbkg = True
#useddbkg = False
labels = {
    "nleps"       : "Nleps",
    "ht"          : "$H_{T}$",
    "met"         : "$p_{T}^{miss}$",
    "mll"         : "$m_{ll}$",
    "njets"       : "Njets",
    "nbtags"      : "Nbtags",    
    "nele"        : "No of electrons", 
    "pt1"         : "$p_T$(lep1)",
    "pt2"         : "$p_T$(lep2)",
    "pt3"         : "$p_T$(lep3)",    
    "eta1"        : r"$\eta$(lep1)",
    "eta2"        : r"$\eta$(lep2)",
    "ptj1"        : "$p_T$ - jet 1",
    "ptj2"        : "$p_T$ - jet 2",
    "ptj3"        : "$p_T$ - jet 3",

    "fwd_jetpt"   : "$p_T$ - fwd jet",    
    "mtmin"       : "$m_{T}^\\mathrm{min}$",
    "mt1"         : "$m_{T}^1$",
    "mt2"         : "$m_{T}^2$",
    "nvtx"        : "# good vertices",    
    "dphil1l2"    : r"$\Delta\phi(l_1,l_2)$",
    "dphil1met"   : r"$\Delta\phi(l_1,p_{T}^{miss})$",

    "dphimetj1"   : r"$\Delta\phi(j_1,p_{T}^{miss})$",
    "drl1l2"      : r"$\Delta R(l_1,l_2)$",
    
    "ptrel1"      : "$p_T^{rel}$(lep1)",
    "ptrel2"      : "$p_T^{rel}$(lep2)",
    "ptratio1"    : "$p_T^{ratio}$(lep1)",
    "ptratio2"    : "$p_T^{ratio}$(lep2)",
    "miniiso1"    : "miniiso1",
    "miniiso2"    : "miniiso2",
    "mindrl1j"    : r"$\Delta R^{min}(l_1,j)$",

    "mindrl1bt"   : r"$\Delta R^{min}(l_1,btag)$",

    "sr_fcnc"     : "Signal Regions",
    
    "l1dxy"       : "$l^{1}_{dxy}$",
    "l2dxy"       : "$l^{2}_{dxy}$",
    "l1dz"        : "$l^{1}_{dz}$",
    "l2dz"        : "$l^{2}_{dz}$",
    "mossf"       : "MOSSF",            



#    "bdt"         :"Event Discriminator",
#    "bdt_hut"     :"Event Discriminator hut",
#    "bdt_hct"     :"Event Discriminator hct",
#    
#    "bdt_hut_ttbar" :"Event Discriminator hut ttbar",
#    "bdt_hut_ttv"   :"Event Discriminator hut ttv",
#    "bdt_hct_ttbar" :"Event Discriminator hct ttbar",
#    "bdt_hct_ttv"   :"Event Discriminator hct ttv",

    #"htb"        : r"$H_{T}$(b-jets)",
    #"nlb40"      : r"N-loose b-tags, $p_{T}>40$",
    #"ntb40"      : r"N-tight b-tags, $p_{T}>40$",
    
}

sig1 = "truthfakes_ttsl"

#sig2 = "fakes_ttsl"
#sig2 = "fakes_mc"

d_label_colors = {
            "dy":                      (r"DY+jets",             [0.4, 0.6, 1.0]),
            "fakes":                   (r"Nonprompt lep.",      [0.85, 0.85, 0.85]),
            "fakes_exp":               (r"Prediction",          [0.85, 0.85, 0.85]),
            "fakes_ttsl":              (r"Prediction",          [0.85, 0.85, 0.85]),
            "fakes_wjets":             (r"Prediction",          [0.85, 0.85, 0.85]),
            "fakes_mc":                (r"Nonprompt lep.",      [0.85, 0.85, 0.85]),
            "fakes_mc_ml":             (r"Nonprompt lep.",      [0.85, 0.85, 0.85]),
            "ttsl":                    (r"ttsl",                [0.85, 0.85, 0.85]),
            "flips":                   (r"Charge misid.",       [0.4, 0.4, 0.4]),
            "flips_mc":                (r"Charge misid.",       [0.4, 0.4, 0.4]),
            "rares":                   ("Rare",                 [1.0, 0.4, 1.0]),
            "singletop":               ("Single Top",           [1.0, 0.4, 0.0]),
            "tt":                      (r"$t\bar{t}$",          [0.8, 0.8, 0.8]),
            "ttfake":                  (r"$t\bar{t}$ Nonprompt",[0.85, 0.85, 0.85]),
            "wjets":                   (r"W+jets",              [113./255,151./255,44./255]),
            "tth":                     (r"$t\bar{t}H$",         [0.4, 0.4, 0.6]),
            "ttw":                     (r"$t\bar{t}W$",         [0.0, 0.4, 0.0]),
            "ttz":                     (r"$t\bar{t}Z$",         [0.4, 0.8, 0.4]),
            #"zbb":                     (r"$t\bar{t}Z$(bb)",     [0.3, 0.6, 0.4]),
            "wz":                      (r"WZ" ,                 [1.0,0.8,0.0]),
            "vv":                      (r"VV",                  [0.0, 0.4, 0.8]),
            "ttvv":                    (r"$t\bar{t}$VV",        [0.0, 0.4, 0.8]),
            "tttt":                    (r"$t\bar{t}t\bar{t}$",  [0.786,0.147,0.022]),
            "wgamma":                  (r"W+$\gamma$",          "#9D7ABF"),
            "zgamma":                  (r"Z+$\gamma$",          "#8154AD"),
            "xg":                      (r"X+$\gamma$",          "#54267F"),
            "raresnoxg":               ("Rare",                 [1.0, 0.4, 1.0]),
        }

d_flat_systematics = {
    "Nonprompt lep."        : 0.4,        
    "Prediction"            : 0.4,        
    "Charge misid."         : 0.3,
    "Rare"                  : 0.5,
    "$t\bar{t}W$"           : 0.3,
    "$t\bar{t}Z$"           : 0.3,
    "$t\bar{t}H$"           : 0.3,
    "X+$\gamma$"            : 0.5,
    "$t\bar{t}W$VV"         : 0.5,
    }

bginfo = {        
    "ssbr"          : { k:d_label_colors[k] for k in [ "fakes_ttsl" ] },
    "ssbrfakeel"    : { k:d_label_colors[k] for k in [ "fakes_ttsl" ] },
    "ssbrfakemu"    : { k:d_label_colors[k] for k in [ "fakes_ttsl" ] },    
}
 
# make these global for multiprocessing since uproot file objects can't be pickled
files, other_files = {}, {}
#print bginfo
def worker(info):
    global files, other_files

    outputdir, year, lumi, region, flav, var = info
    title = region.upper()
    xlabel = labels[var]
    hname = "{}_{}_{}".format(region,var,flav)

    if other_files:
        bgs = [
            sum([Hist1D(files[proc][hname],label=label,color=color)] + [Hist1D(other_files[y][proc][hname],label=label,color=color) for y in other_files.keys()])
            for proc,(label,color) in sorted(bginfo[region].items())
            ]
        data = sum([Hist1D(files["data"][hname])] + [Hist1D(other_files[y]["data"][hname]) for y in other_files.keys()])
        
        sigs = [
            sum([Hist1D(files[sig1][hname],label="Fakes", color = "#9D7ABF" )] + [Hist1D(other_files[y][sig1][hname],label="Fakes", color = "#9D7ABF") for y in other_files.keys()]),
            
            ]        
        
    else:
        bgs = [Hist1D(files[proc][hname], label=label,color=color) for proc,(label,color) in sorted(bginfo[region].items())]
        data = Hist1D(files["data"][hname])
        sigs = [Hist1D(files[sig1][hname],label="Fakes", color = "#9D7ABF") ]           
        

    data.set_attr("label", "Data [{}]".format(int(data.get_integral())))
    sigs[0].set_attr("label", "Fakes [{:.1f}]".format(sigs[0].get_integral()))
    #sigs[1].set_attr("label", "hct [{:.1f}]".format(sigs[1].get_integral()))
    #sum(sigs).set_attr("color", [1.0, 0.4, 1.0])
    
    #print "S1"
    #print "data ",data.get_integral()
    #print "bgs ", abs(sum(bgs).get_integral())
    #if data.get_integral() < 1e-6: return
    if abs(sum(bgs).get_integral()) < 1e-6: return

    #print "S2"
    do_bkg_syst = True
    
    bgs = sorted(bgs, key=lambda bg: bg.get_integral())
    sf = data.get_integral()/sum(bgs).get_integral()
    #bgs = [bg*sf for bg in bgs]
    # bgs = [bg*1 for bg in bgs]
    

    for bg in bgs:
        # add flat systematic to stat unc in quadrature                      
        #print bg.get_attr("label")
        bg._errors = np.hypot(bg._counts*d_flat_systematics.get(bg.get_attr("label"),0.),bg._errors)

    if plotdata:
        title += " data/MC={:.2f}".format(sf)

    if other_files:
        fname = "{}/run2_{}_{}_{}.pdf".format(outputdir,region,var,flav)
    else:
        fname = "{}/year{}_{}_{}_{}.pdf".format(outputdir,year,region,var,flav)
    
    #print "S3"
    if plotdata :
        plot_stack(bgs=bgs, 
               data=data, 
               sigs = sigs,
               title=title,                
               xlabel=xlabel, 
               filename=fname,
               cms_type = "Preliminary",
               # do_log=True,
               do_bkg_syst=do_bkg_syst,
               lumi = lumi,
               ratio_range=[0.0,2.0],
               mpl_title_params=dict(fontsize=(8 if len(str(lumi))>=5 else 9)),
               # ratio_range=[0.5,1.5],
               )

        fname_log = fname.replace(".pdf","_log.pdf").replace(".png","_log.png")
        plot_stack(bgs=bgs, 
               data=data, 
               sigs = sigs, 
               title=title,                
               xlabel=xlabel, 
               filename=fname_log,
               cms_type = "Preliminary",
               do_log=True,
               do_bkg_syst=do_bkg_syst,
               lumi = lumi,
               ratio_range=[0.0,2.0],
               mpl_title_params=dict(fontsize=(8 if len(str(lumi))>=5 else 9)),
               # ratio_range=[0.5,1.5],
               )
    if not plotdata:
        plot_stack(bgs=bgs, 
                   #data=data, 
                   sigs = sigs,
                   ratio = sigs[0].divide(sum(bgs)),
                   title=title, 
                   xlabel=xlabel, 
                   filename=fname,
                   cms_type = "Preliminary",
                   # do_log=True,
                   do_bkg_syst=do_bkg_syst,
                   lumi = lumi,
                   ratio_range=[0.0,2.0],
                   mpl_title_params=dict(fontsize=(8 if len(str(lumi))>=5 else 9)),
                   mpl_ratio_params={"label":"hut/Bkgd"},
                   # ratio_range=[0.5,1.5],
                   )

        fname_log = fname.replace(".pdf","_log.pdf").replace(".png","_log.png")
        plot_stack(bgs=bgs, 
               #data=data, 
               sigs = sigs, 
               ratio = sigs[0].divide(sum(bgs)),
               title=title, 
               xlabel=xlabel, 
               filename=fname_log,
               cms_type = "Preliminary",
               do_log=True,
               do_bkg_syst=do_bkg_syst,
               lumi = lumi,
               ratio_range=[0.0,2.0],
               mpl_title_params=dict(fontsize=(8 if len(str(lumi))>=5 else 9)),
               mpl_ratio_params={"label":"hut/Bkgd"},
               # ratio_range=[0.5,1.5],
               )
    # os.system("ic {}".format(fname))
    #write_table(data,bgs,signal=sigs,show_errors=False,outname=fname.replace(".pdf",".txt"))
    return fname

def make_plots(outputdir="plots", inputdir="outputs", year=2016, lumi="35.9", other_years=[], regions=[], flavs=["ee","em","mm","in"]):
    global files, other_files

    os.system("mkdir -p {}/".format(outputdir))

    files = { proc:uproot.open("{}/histos_{}_{}.root".format(inputdir,proc,year)) for proc in (list(set(sum(map(lambda x:x.keys(),bginfo.values()),[])))+ ["data"] + [sig1]) }
    other_files = {}
    for y in other_years:
        other_files[y] = { proc:uproot.open("{}/histos_{}_{}.root".format(inputdir,proc,y)) for proc in (list(set(sum(map(lambda x:x.keys(),bginfo.values()),[]))) + ["data"] + [sig1]) }
        
    #print files
    # for region in ["htnb1mc","htnb1","os","osloose","br","crw","crz","tt_isr_reweight_check"]:
    # regions = ["htnb1mc","htnb1","htnb1mcmu","htnb1mu","os","os_noht","osloose","br","crw","crz"]
    regions = regions or ["ssbr","mlbr"]
    flavs = flavs or ["ee","em","mm","in"]
    varss = labels.keys()
    infos = [[outputdir,year,lumi]+list(x) for x in itertools.product(regions,flavs,varss)]

    os.nice(10)
    from multiprocessing import Pool as ThreadPool
    pool = ThreadPool(25)
    for res in pool.imap_unordered(worker,infos):
        if res:
            print "Wrote {}".format(res)


if __name__ == "__main__":
    ## 
    #if useddbkg:        
    regions = [
        
        "ssbr",
        "ssbrfakeel",
        "ssbrfakemu"
        
    ]
               
    #flavs = ["in","ee","em","mm"]
    flavs = ["in"]
    
    #inputdir = "outputs_v3p31_apr27/"
    inputdir = "outputs_v3p31_may6/"
    outputdir = inputdir+"plots"

    #print regions

    make_plots(
        outputdir=outputdir,
        inputdir=inputdir,
        regions = regions, flavs = flavs,
        year=2016,
        lumi="35.9",
    )
    
    # 2017 alone
    make_plots(
        outputdir=outputdir,
        inputdir=inputdir,
        regions = regions, flavs = flavs,
        year=2017,
        lumi="41.5",
    )

    # 2018 alone
    make_plots(
        outputdir=outputdir,
        inputdir=inputdir,
        regions = regions, flavs = flavs,
        year=2018,
        lumi="59.7",
    )

    # # 2016 + 2018 + 2017
    make_plots(
        outputdir=outputdir,
        inputdir=inputdir,
        regions = regions, flavs = flavs,
        year=2016, 
        lumi="137",
        other_years = [2017,2018],
    )
    
    #os.system("outputs_v3p31_SR_FCNC_Jan8/sr_plot v3p31_sr_plot_Jan8 ")
