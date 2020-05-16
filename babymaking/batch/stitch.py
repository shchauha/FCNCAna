import os
import sys
import ROOT as r

import monkeyroot

# f1 = r.TFile("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_fakesv2/output/year_2016/TTBAR_PH.root")
# f2 = r.TFile("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_fakesv2/output/year_2016/TTDLht500.root")
f1 = r.TFile("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.10_prefire/output/year_2016/TTBAR_PH.root")
# f2 = r.TFile("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.10_prefire/output/year_2016/TTDLht500.root")
f2 = r.TFile("output_ttslht500.root")
t1 = f1.Get("t")
t2 = f2.Get("t")

# t3 = r.TChain("t")
# t3.Add("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_lowpt3/output/year_2016/TTWnlo.root")
# t3.Add("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_lowpt3/output/year_2016/TTZnlo.root")
# t3.Add("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_lowpt3/output/year_2016/TTHtoNonBB.root")

# t4 = r.TChain("t")
# t4.Add("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_fakesv2/output/year_2016/TTDL.root")
# t4.Add("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_lowpt3/output/year_2016/TTSLtop.root")
# t4.Add("/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.09_lowpt3/output/year_2016/TTSLtopbar.root")

c1 = r.TCanvas()

# NOTE DL
# h1 = r.TH1F("h1","h1",15,0,15)
# h2 = r.TH1F("h2","h2",15,0,15)
# t1.Draw("ngenjets30>>h1","scale1fb*(genht30>500 && nleptonicW==2)","goff")
# t2.Draw("ngenjets30>>h2","scale1fb*(genht30>500 && nleptonicW==2)","goff")
# sf = h1.Integral(h1.FindBin(7),-1)/h2.Integral(h2.FindBin(7),-1)
# print "sf =",sf
# h2.Scale(sf)
# ve1 = h1.GetBinValueErrors()
# ve2 = h2.GetBinValueErrors()
# print map(lambda x: (x[1]/x[0]).round(3), zip(ve1[8:],ve2[8:]))

# NOTE SL
h1 = r.TH1F("h1","h1",15,0,15)
h2 = r.TH1F("h2","h2",15,0,15)
t1.Draw("ngenjets30>>h1","scale1fb*(genht30>500 && nleptonicW==1)","goff")
t2.Draw("ngenjets30>>h2","scale1fb*(genht30>500 && nleptonicW==1)","goff")
sf = h1.Integral(h1.FindBin(9),-1)/h2.Integral(h2.FindBin(9),-1)
print "sf =",sf # 
h2.Scale(sf)
ve1 = h1.GetBinValueErrors()
ve2 = h2.GetBinValueErrors()
print ve1[10:], ve2[10:]
print map(lambda x: (x[1]/x[0]).round(3), zip(ve1[10:-1],ve2[10:-1]))

# h1 = r.TH1F("h1","h1",15,0,15)
# h2 = r.TH1F("h2","h2",15,0,15)
# # h1 = r.TH2F("h1","h1",2,-0.5,1.5,2,-0.5,1.5)
# # h2 = r.TH2F("h2","h2",2,-0.5,1.5,2,-0.5,1.5)
# # h1.SetMarkerSize(2)
# # h2.SetMarkerSize(2)
# t1.Draw("ngenjets30>>h1","scale1fb*(passfilter && genht30>500 && nleptonicW==2)","goff")
# t2.Draw("ngenjets30>>h2","scale1fb*(passfilter && genht30>500 && nleptonicW==2)","goff")
# # t1.Draw("ngenjets30>=7:genht30>500>>h1","scale1fb*(nleptonicW==2)","goff")
# # t2.Draw("ngenjets30>=7:genht30>500>>h2","scale1fb*(nleptonicW==2)","goff")
# # t2.Draw("ngenjets30>=7:genht30>500>>h2","scale1fb*(passfilter && nleptonicW==2)","goff")

# h1 = r.TH1F("h1","h1",50,0,1000)
# h2 = r.TH1F("h2","h2",50,0,1000)
# # t1.Draw("genht30>>h1","scale1fb*(passfilter && ngenjets30>=7 && nleptonicW==2)","goff")
# # t2.Draw("genht30>>h2","scale1fb*(passfilter && ngenjets30>=7 && nleptonicW==2)","goff")
# t1.Draw("genht30>>h1","scale1fb*(nleptonicW==2 && passfilter)","goff")
# t2.Draw("genht30>>h2","scale1fb*(nleptonicW==2 && passfilter)","goff")

# # h1 = r.TH1F("h1","h1",15,0,15)
# # h2 = r.TH1F("h2","h2",15,0,15)
# h1 = r.TH1F("h1","h1",50,0,1000)
# h2 = r.TH1F("h2","h2",50,0,1000)
# t1.Draw("genht30>>h1","scale1fb*(passfilter && ngenjets30>=7 && nleptonicW==2)","goff")
# t2.Draw("genht30>>h2","scale1fb*(passfilter && ngenjets30>=7 && nleptonicW==2)","goff")
# # t1.Draw("genht30>>h1","scale1fb*(ngenjets>=7 && genht30>=500 && nleptonicW==2)","goff")
# # t2.Draw("genht30>>h2","scale1fb*(ngenjets>=7 && genht30>=500 && nleptonicW==2)","goff")
# h1.Scale(1./h1.Integral())
# h2.Scale(1./h2.Integral())



# don't take inclusive sample if 
# ngenjets30>=7 && nleptonicW==2 && genht30>500
# and take ht500 sample instead

# h1 = r.TH1F("h1","h1",16,0.5,16.5)
# h2 = r.TH1F("h2","h2",16,0.5,16.5)
# # h3 = r.TH1F("h3","h3",16,0.5,16.5)
# # h4 = r.TH1F("h4","h4",16,0.5,16.5)
# t1.SetAlias("keep", "!(ngenjets30>=7 && nleptonicW==2 && genht30>500)")
# t2.SetAlias("keep", " (ngenjets30>=7 && nleptonicW==2 && genht30>500)")
# # t1.SetAlias("keep", "passfilter")
# # t2.SetAlias("keep", "passfilter")
# # for t in [t1,t2,t3,t4]:
# for t in [t1,t2]:
#     t.SetAlias("hasfake", "(lep1_motherID<=0) + (lep2_motherID<=0) + (nleps>2)*(lep3_motherID<=0)")
# print "EHRE"
# # t1.Draw("sr-2>>h1","35.9*scale1fb*{}*(keep && br && hyp_class==3 && hasfake>0 && sr>1)".format(1.0),"goff")
# # t2.Draw("sr-2>>h2","35.9*scale1fb*{}*(keep && br && hyp_class==3 && hasfake>0 && sr>1)".format(1.0),"goff")
# t1.Draw("sr-2>>h1","35.9*scale1fb*{}*(keep && br && hyp_class==3 && sr>1)".format(1.0),"goff")
# t2.Draw("sr-2>>h2","35.9*scale1fb*{}*(keep && br && hyp_class==3 && sr>1)".format(1.0),"goff")
# # t3.Draw("sr-2>>h3","35.9*scale1fb*{}*(br && hyp_class==3 && hasfake>0 && sr>2)".format(1.0),"goff")
# # t4.Draw("sr-2>>h4","35.9*scale1fb*{}*(br && hyp_class==3 && hasfake>0 && sr>2)".format(1.0),"goff")
# h1.Add(h2)
# # h1.Add(h3)
# # print list(h3)
# h1.SetMinimum(0.001)
# h2.SetMinimum(0.001)
# c1.SetLogy()

h1.MakePretty(r.kAzure+2)
h2.MakePretty(r.kRed)
# h3.MakePretty(r.kOrange+2)
h1.Draw("histe")
# h1.SetMinimum(0.)
# h4.Draw("histesame")
h2.Draw("samehiste")
# h3.Draw("samehiste")
c1.SaveAs("test.pdf")
os.system("ic test.pdf")

# # h1.MakePretty(r.kAzure+2)
# # h2.MakePretty(r.kRed)
# # h3.MakePretty(r.kOrange+2)
# # h1.Draw("colztext")
# h2.Draw("colztext")
# # h1.SetMinimum(0.)
# # h4.Draw("histesame")
# # h2.Draw("samehiste")
# # h3.Draw("samehiste")
# c1.SaveAs("test.pdf")
# os.system("ic test.pdf")
