import ROOT as r
# r.gROOT.ProcessLine(".L /home/users/namin/2018/fourtop/all/FTAnalysis/analysis/misc/year_2018/isr_weights.h ") # FIXME

# 0. Derive isr weights from 0/1 extra parton samples in $FTBASE/analysis/crplots
# 1. Put appropriate isr-reweighting function in process line
# 2. Comment out any scales used in the function since we are going to derive them here
# 3. Edit basedir and sname/snum below
# 4. Run the script and put the sfs into the header file
# 5. Now, re-running this script should give you 1.0 in principle

# r.gROOT.ProcessLine(".L /home/users/namin/2018/fourtop/all/FTAnalysis/analysis/misc/year_2017/isr_weights.h ")
# basedir = "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.21_fix2017/output/year_2017/" # FIXME

r.gROOT.ProcessLine(".L /home/users/namin/2018/fourtop/all/FTAnalysis/analysis/misc/year_2018/isr_weights.h ")
basedir = "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.23/output/year_2018/" # FIXME

# r.gROOT.ProcessLine(".L /home/users/namin/2018/fourtop/all/FTAnalysis/analysis/misc/year_2016/isr_weights.h ")
# basedir = "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.23/output/year_2016/"

for sname,snum in [
        ("TTWnlo",1),
        ("TTZnlo",2),
        ("TTBAR*",10),
        ]:
    ch = r.TChain("t")
    ch.Add("{}/{}.root".format(basedir,sname))
    extra = ""
    hw = r.TH1D("hw","",1,0,1)
    hunw = r.TH1D("hunw","",1,0,1)
    ch.Draw("0.5>>hw","isrWeight(nisrMatch, {}) * scale1fb {}".format(snum,extra))
    ch.Draw("0.5>>hunw","1 * scale1fb {}".format(extra))
    hunw = r.gDirectory.Get("hunw")
    hw = r.gDirectory.Get("hw")
    vhunw = hunw.Integral()
    vhw = hw.Integral()
    sf = vhunw / vhw
    print sname, snum, sf
