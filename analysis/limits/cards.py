#!/usr/bin/env python
from ROOT import *
from array import *
gROOT.SetBatch(True)
import sys, os
import math
from fractions import Fraction
#ss_nbtags[1]

def BinContents( hist,bin):
    nbinx  = hist.GetNbinsX()
    #print "nbinx ", nbinx    
    if(bin<nbinx):
        content = hist.GetBinContent(bin)
        if(content<0.0001):
            content = 0.0001
    if(bin==hist.GetNbinsX()):
        content =  hist.GetBinContent(bin)+hist.GetBinContent(bin+1)
        if(content<0.0001):
            content = 0.0001
    return content
            
def BinError(hist, bin):
    nbinx  = hist.GetNbinsX()
    #print "nbinx ", nbinx
    if(bin<nbinx):
        error = hist.GetBinError(bin)
        if(error<0.0001):
            error = 0.0001
    
    if(bin==nbinx):
        error =pow(hist.GetBinError(bin)*hist.GetBinError(bin)+hist.GetBinError(bin+1)*hist.GetBinError(bin+1),0.5)
        if(error<0.0001):
            error = 0.0001        
    return error

outdir = "CNC"
os.system("mkdir -p {}".format(outdir))
def createDatcards(varname="Variable", signalname="Signal", srname="SR", year="2016", debug=False) :
    #basedir = "/home/shchauha/2019/Analysis/FCNCAna/analysis/yields/v3.24/"
    #basedir = "/home/users/shchauha/2019/FCNCAna/analysis/yields/v3.24_test/"
    bginput = [basedir+"histos_fcnc_"+year+".root",
               basedir+"histos_ttdl_"+year+".root",
               basedir+"histos_ttsl_"+year+".root",
               basedir+"histos_tth_"+year+".root",
               basedir+"histos_ttw_"+year+".root",
               basedir+"histos_ttz_"+year+".root",
               basedir+"histos_wz_"+year+".root",
               #basedir+"histos_dy.root"        
    ]
    sginput = [basedir+"histos_fcnc_"+year+".root"]
    datainput = [basedir+"histos_fcnc_"+year+".root"]

    fd = TFile(basedir+"histos_fcnc_"+year+".root","read")
    fs = TFile(basedir+"histos_fcnc_"+year+".root","read")
    hdata = fd.Get(varname).Clone()
    hs = fs.Get(varname).Clone()
    
    #print "signal events ",hs.GetBinContent(3), " my func ", BinContents(hs,3)
    #print "signal events err ",hs.GetBinError(3), " my func ", BinError(hs,3)
    
    process = [signalname,"chargeflip","nonprompt","tth","ttw","ttz","wz"]
    #process.append(signalname)
    bgname =[] 
    #hbg = []

    flatsystematics=[
        ("WZNORM"    ,0.30,               "wz"),
        ("TTWNORM"   ,0.30,               "ttw"),
        ("TTZNORM"   ,0.30,               "ttz"),
        ("TTHNORM"   ,0.30,               "tth"),
        ("nonprompt" , 0.30,              "nonprompt"),
        ("chargeflip", 0.30,              "chargeflip"),
        ("IDISO"     , 0.06,              signalname,"ttw","ttz","tth"), 
        ("LUMI"      ,0.024,              signalname,"ttw","ttz","tth"),
    ]
    
    for i in range(len(bginput)):
        bgname.append(bginput[i].replace("_"+year+".root","").replace(basedir+"histos_",""))
    
    for bin in range(hs.GetNbinsX()+1):
        if(hs.GetBinContent(bin)==0):
            print "skipped bin ", bin
            continue
        #print "bin content ", hs.GetBinContent(bin)
        outname = outdir+"/"+year+"_"+srname+"_"+str(bin)+".txt"
        ustr = srname+"_"+year+"_"+str(bin)
        f = open(outname, 'w')
        orig_stdout = sys.stdout
        sys.stdout = f
    
        print "imax ",1
        print "jmax ",len(bginput)-1
        print "kmax ",13

        print "bin "+ustr
        print "observation",    0

        print "#----------------"
        print "bin     ",(ustr+" ")*(len(bginput))
        print "process",
        for i in range(len(bgname)):
            print i,
        print ""
        print "process ",
        for i in range(len(bgname)):
            print bgname[i],
        print ""
        print "rate    ", 
        for i in range(len(bgname)):            
            temp = TFile(bginput[i],"read")
            #print "varname ",varname
            temph = temp.Get(varname).Clone()
            print BinContents(temph,bin),

        print ""    
        print "#----------------"
        for i in range(len(bgname)):            
            print bgname[i]+"_"+ustr+"_STAT gmN ",
            #print bgname[i]+"_STAT gmN ",
            temp = TFile(bginput[i],"read")
            temph = temp.Get(varname).Clone()
            bincontent = BinContents(temph,bin)
            binerror = BinError(temph,bin)
            relbinerror = binerror/bincontent
            processn=int(math.ceil(pow(1./relbinerror,2)))
            processweight=bincontent/float(processn)
            print processn,"  -  "*(i), processweight,"  -  "*(len(bgname)-i-1)


        print "ttdl_SYST lnN   -   1.30  -    -    -    -    -   "
        print "ttsl_SYST lnN   -    -   1.40  -    -    -    -  "    
        print "tth_SYST lnN   -   -  - 1.30    -    -    -     "
        print "ttw_SYST lnN   -   -    -  -  1.30  -    -   "
        print "ttz_SYST lnN   -   -    -  - -   1.30     -  "
        print "wz_SYST lnN   -    -    -  -  - - 1.30    "
        #print "DY_high_SYST lnN   -    -  -  -    -    - -   1.20  "    
    
        sys.stdout = orig_stdout
        f.close()        
        
if __name__ == '__main__':
    
    basedir          = sys.argv[1]
    varname          = sys.argv[2]
    signalname       = sys.argv[3]
    srname           = sys.argv[4]
    year             = sys.argv[5]
    #debugoff         = int(sys.argv[5])
    
    #debug=True
    
    #if debugoff==1 : debug=False
    debug=False

    if debug : print "\npython dataCardProducer.py basedir Variable  Signal SRname year debugOff \n"
    createDatcards(varname, signalname, srname, year, debug)
    exit(0)
