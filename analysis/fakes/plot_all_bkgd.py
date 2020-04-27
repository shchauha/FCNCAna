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
    "ht"          : "$H_{T}$",
    "met"         : "$p_{T}^{miss}$",
    "mll"         : "$m_{ll}$",
    "njets"       : "Njets",
    "nbtags"      : "Nbtags",
    "nbtagszoom"  : "Nbtags",
    "nele"        : "No of electrons", 
    "pt1"         : "$p_T$(lep1)",
    "pt2"         : "$p_T$(lep2)",
    "pt3"         : "$p_T$(lep3)",    
    "eta1"        : r"$\eta$(lep1)",
    "eta2"        : r"$\eta$(lep2)",
    "ptj1"        : "$p_T$ - jet 1",
    "ptj2"        : "$p_T$ - jet 2",
    "ptj3"        : "$p_T$ - jet 3",
    "ptbt1"       : "$p_T$ - btag 1",
    "ptbt2"       : "$p_T$ - btag 2",
    "fwd_jetpt"   : "$p_T$ - fwd jet",    
    "mtmin"       : "$m_{T}^\\mathrm{min}$",
    "mt1"         : "$m_{T}^1$",
    "mt2"         : "$m_{T}^2$",
    "nvtx"        : "# good vertices",    
    "dphil1l2"    : r"$\Delta\phi(l_1,l_2)$",
    "dphil1met"   : r"$\Delta\phi(l_1,p_{T}^{miss})$",
    "dphil2met"   : r"$\Delta\phi(l_2,p_{T}^{miss})$",
    "dphimetj1"   : r"$\Delta\phi(j_1,p_{T}^{miss})$",
    "drl1l2"      : r"$\Delta R(l_1,l_2)$",
    
    "ptrel1"      : "$p_T^{rel}$(lep1)",
    "ptrel2"      : "$p_T^{rel}$(lep2)",
    "ptratio1"    : "$p_T^{ratio}$(lep1)",
    "ptratio2"    : "$p_T^{ratio}$(lep2)",
    "miniiso1"    : "miniiso1",
    "miniiso2"    : "miniiso2",
    "mindrl1j"    : r"$\Delta R^{min}(l_1,j)$",
    "mindrl2j"    : r"$\Delta R^{min}(l_2,j)$",
    "mindrl1bt"   : r"$\Delta R^{min}(l_1,btag)$",
    "mindrl2bt"   : r"$\Delta R^{min}(l_2,btag)$",
    
    "l1dxy"       : "$l^{1}_{dxy}$",
    "l2dxy"       : "$l^{2}_{dxy}$",
    "l1dz"        : "$l^{1}_{dz}$",
    "l2dz"        : "$l^{2}_{dz}$",
    "mossf"       : "MOSSF",            
    "bdt"         :"Event Discriminator",
    "bdt_hut"     :"Event Discriminator hut",
    "bdt_hct"     :"Event Discriminator hct",

    #"htb"        : r"$H_{T}$(b-jets)",
    #"nlb40"      : r"N-loose b-tags, $p_{T}>40$",
    #"ntb40"      : r"N-tight b-tags, $p_{T}>40$",
    
}

d_label_colors = {
            "dy":                      (r"DY+jets",             [0.4, 0.6, 1.0]),
            "fakes":                   (r"Nonprompt lep.",      [0.85, 0.85, 0.85]),
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
    "Charge misid."         : 0.3,
    "Rare"                  : 0.5,
    "$t\bar{t}W$"           : 0.3,
    "$t\bar{t}Z$"           : 0.3,
    "$t\bar{t}H$"           : 0.3,
    "X+$\gamma$"            : 0.5,
    "$t\bar{t}W$VV"         : 0.5,
    }

bginfodata = {}
bginfomc = {}
bginfoall = {}

bginfoall= {        
    "ssbr"          : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "ss0b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "ss1b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "ss2b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "mlbr"          : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "fakes_mc_ml", "xg", "ttvv", "rares" ] },                    
    "ml1b1j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "fakes_mc_ml", "xg", "ttvv", "rares" ] },                    
    "ml2b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "fakes_mc_ml", "xg", "ttvv", "rares" ] },                    
    "mlbrinc"       : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "fakes_mc_ml", "xg", "ttvv", "rares" ] },                    
    
    "mllowmetonz2b" : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "fakes_mc_ml", "xg", "ttvv", "rares" ] },                    
    }



bginfodata = {        
    "ssbr"          : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "xg", "ttvv", "rares"] },
    "ss0b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "xg", "ttvv", "rares"] },
    "ss1b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "xg", "ttvv", "rares"] },
    "ss2b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "flips", "xg", "ttvv", "rares"] },
    "mlbr"          : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "xg", "ttvv", "rares" ] },                    
    "ml1b1j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "xg", "ttvv", "rares" ] },                    
    "ml2b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "xg", "ttvv", "rares" ] },                    
    "mlbrinc"       : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "xg", "ttvv", "rares" ] },                    
    
    "mllowmetonz2b" : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "xg", "ttvv", "rares" ] },                    
    }

bginfomc = {
    "ssbr"          : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "ss0b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "ss1b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "ss2b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc", "flips_mc", "xg", "ttvv", "rares"] },
    "mlbr"          : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc_ml", "xg", "ttvv", "rares" ] },        
    "ml1b1j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc_ml", "xg", "ttvv", "rares" ] },            
    "ml2b2j"        : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc_ml", "xg", "ttvv", "rares" ] },            
    "mlbrinc"       : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc_ml", "xg", "ttvv", "rares" ] },            
    "osbr"          : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc", "flips_mc", "xg", "ttvv", "rares" ] },
    "tl"            : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes_mc", "flips_mc", "xg", "ttvv", "rares" ] },      
    "mllowmetonz2b" : { k:d_label_colors[k] for k in [ "ttw", "tth", "ttz", "fakes", "xg", "ttvv", "rares" ] },                    
    }    
# make these global for multiprocessing since uproot file objects can't be pickled
files, other_files = {}, {}

def worker(info):
    global files, other_files

    outputdir, year, lumi, region, flav, var = info
    title = region.upper()
    xlabel = labels[var]
    hname = "{}_{}_{}".format(region,var,flav)

    if other_files:
        bgs = [
            sum([Hist1D(files[proc][hname],label=label,color=color)] + [Hist1D(other_files[y][proc][hname],label=label,color=color) for y in other_files.keys()])
            for proc,(label,color) in sorted(bginfodata[region].items())
            ]
        data = sum([Hist1D(files["data"][hname])] + [Hist1D(other_files[y]["data"][hname]) for y in other_files.keys()])
        sigs = [
            sum([Hist1D(files[proc][hname],label=label,color=color)] + [Hist1D(other_files[y][proc][hname],label=label,color=color) for y in other_files.keys()])
            for proc,(label,color) in sorted(bginfomc[region].items())
            ]

#        sigs = [
#            sum([Hist1D(files["fcnc_hut"][hname],label="hut", color = "#9D7ABF" )] + [Hist1D(other_files[y]["fcnc_hut"][hname],label="hut", color = "#9D7ABF") for y in other_files.keys()]),
#            sum([Hist1D(files["fcnc_hct"][hname],label="hct", color = "#8154AD" )] + [Hist1D(other_files[y]["fcnc_hct"][hname],label="hct", color = "#8154AD") for y in other_files.keys()])
#            ]        
        

    else:
        bgs = [Hist1D(files[proc][hname], label=label,color=color) for proc,(label,color) in sorted(bginfodata[region].items())]
        data = Hist1D(files["data"][hname])
        sigs = [Hist1D(files[proc][hname], label=label,color=color) for proc,(label,color) in sorted(bginfomc[region].items())]
        #sigs = [Hist1D(files["fcnc_hut"][hname],label="hut", color = "#9D7ABF"), Hist1D(files["fcnc_hct"][hname],label="hct", color = "#8154AD" ) ]           
        
    #print sigs
    data.set_attr("label", "Data [{}]".format(int(data.get_integral())))
    #    sigs[1].set_attr("label", "hct [{:.1f}]".format(sigs[1].get_integral()))
    num = []
    den = []
    for bg in bgs:
        if bg.get_attr("label") == "Nonprompt lep." or bg.get_attr("label") == "Charge misid.":
            num.append(bg)

    for sig in sigs:
        #print sig.get_attr("label")
        if sig.get_attr("label") == "Nonprompt lep." or sig.get_attr("label") == "Charge misid.":        
            den.append(sig)            
    #sigs = [sum(sigs)]
    sigs[0].set_attr("label", "Bkgd (MC) [{:.1f}]".format(sum(sigs).get_integral()))        
    sigs[0].set_attr("color", "#8154AD")    
    #sum(sigs).set_attr("color", [1.0, 0.4, 1.0])    
    if data.get_integral() < 1e-6: return
    if abs(sum(bgs).get_integral()) < 1e-6: return

    do_bkg_syst = True
    
    bgs = sorted(bgs, key=lambda bg: bg.get_integral())
    sf = sum(bgs).get_integral()/sum(sigs).get_integral()
    sf2 = sum(num).get_integral()/sum(den).get_integral()
    #bgs = [bg*sf for bg in bgs]
    # bgs = [bg*1 for bg in bgs]
    
    for bg in bgs:
        # add flat systematic to stat unc in quadrature                      
        #print bg.get_attr("label")
        bg._errors = np.hypot(bg._counts*d_flat_systematics.get(bg.get_attr("label"),0.),bg._errors)

        
    title += " data/MC={:.2f}".format(sf2)

    if other_files:
        fname = "{}/run2_{}_{}_{}.pdf".format(outputdir,region,var,flav)
    else:
        fname = "{}/year{}_{}_{}_{}.pdf".format(outputdir,year,region,var,flav)
        
    if plotdata :
        plot_stack(bgs=bgs, 
               data=data, 
               sigs = [sum(sigs)],
               ratio = sum(num).divide(sum(den)),                                  
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
               sigs = [sum(sigs)], 
               ratio = sum(num).divide(sum(den)),                                  
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

    #print len(den)    
    if not plotdata:
        plot_stack(bgs=bgs, 
                   #data=data, 
                   sigs = [sum(sigs)],
                   ratio = sum(num).divide(sum(den)),                                  
                   #ratio = (num[0]+num[1]).divide(den[0]+den[1]),
                   title=title, 
                   xlabel=xlabel, 
                   filename=fname,
                   cms_type = "Preliminary",
                   # do_log=True,
                   do_bkg_syst=do_bkg_syst,
                   lumi = lumi,
                   ratio_range=[0.0,2.0],
                   mpl_title_params=dict(fontsize=(8 if len(str(lumi))>=5 else 9)),
                   mpl_ratio_params={"label":"Data/MC"},
                   # ratio_range=[0.5,1.5],
                   )

        fname_log = fname.replace(".pdf","_log.pdf").replace(".png","_log.png")
        plot_stack(bgs=bgs, 
               #data=data, 
               sigs = [sum(sigs)], 
               ratio = sum(num).divide(sum(den)),               
               #ratio = (num[0]+num[1]).divide(den[0]+den[1]),
               title=title, 
               xlabel=xlabel, 
               filename=fname_log,
               cms_type = "Preliminary",
               do_log=True,
               do_bkg_syst=do_bkg_syst,
               lumi = lumi,
               ratio_range=[0.0,2.0],
               mpl_title_params=dict(fontsize=(8 if len(str(lumi))>=5 else 9)),
               mpl_ratio_params={"label":"Data/MC"},
               #mpl_ratio_params={"label":"DDBkgd/MCBkgd"},

               # ratio_range=[0.5,1.5],
               )


    # os.system("ic {}".format(fname))
    #write_table(data,bgs,outname=fname.replace(".pdf",".txt"))
    return fname

def make_plots(outputdir="plots", inputdir="outputs", year=2016, lumi="35.9", other_years=[], regions=[], flavs=["ee","em","mm","in"]):
    global files, other_files

    os.system("mkdir -p {}/".format(outputdir))

    files = { proc:uproot.open("{}/histos_{}_{}.root".format(inputdir,proc,year)) for proc in (list(set(sum(map(lambda x:x.keys(),bginfoall.values()),[])))+["data"]+["fcnc_hut"]+["fcnc_hct"]) }
    other_files = {}
    for y in other_years:
        other_files[y] = { proc:uproot.open("{}/histos_{}_{}.root".format(inputdir,proc,y)) for proc in (list(set(sum(map(lambda x:x.keys(),bginfoall.values()),[])))+["data"]+["fcnc_hut"]+["fcnc_hct"]) }
        
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
    if useddbkg:        
        regions = [
            "ssbr",
            "ss0b2j",
            "ss1b2j",
            "ss2b2j",
            "mlbr",         
            #"mllowmetonz2b",
            "ml2b2j",
            "ml1b1j",
            "mlbrinc",
            #
            ]
        
    if not useddbkg:        
        regions = [
            #"ssbr",
            #"ss0b2j",
            #"ss1b2j",
            #"ss2b2j",
            #"mlbr",         
            #"mllowmetonz2b",
            #"ml2b2j",
            #"ml1b1j",
            #"mlbrinc",
            "osbr",
            "tl",
            "mllowmetonz2b",
            ]
       
    flavs = ["in"]
    # flavs = ["ee","em","mm","in"]
    inputdir = "outputs_v3p31_BDT/"
    outputdir = inputdir+"plots_sr_datamc2"

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
            lumi="137.2",
            other_years = [2017,2018],
            )

    os.system("niceplots outputs_v3p31_BDT/plots_sr_datamc plots_v3.31_sr_datamc")
