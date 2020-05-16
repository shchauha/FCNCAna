import os
import sys
import glob
import analysis.utils.pyrun as pyrun
import ROOT as r
import numpy as np

def embed_ratios(fname):
    f1 = r.TFile(fname,"update")
    hnames = [x.GetName() for x in f1.GetListOfKeys()]
    dens = {}
    for wp in ["med","tight","loose"]:
        nums = {}
        for flav in ["b","c","udsg"]:
            haveeff = len([hn for hn in hnames if ("_{}_Eff_{}".format(wp,flav) in hn)])==1
            if haveeff:
                print "Already have {} efficiency histogram for {} in {}".format(wp,flav,fname)
                continue
            numname = next(hn for hn in hnames if ("_{}_Num_{}".format(wp,flav) in hn))
            nums[flav] = f1.Get(numname)
            if flav not in dens:
                denname = next(hn for hn in hnames if ("_Denom_{}".format(flav) in hn))
                dens[flav] = f1.Get(denname)

        if nums:
            for flav in nums.keys():
                eff = nums[flav].Clone(nums[flav].GetName().replace("_Num","_Eff"))
                eff.Divide(nums[flav],dens[flav],1,1,"B")
                eff.Write()
                print "Wrote {}".format(eff.GetName())
    f1.Close()

def make_two_year_ratio(fname1,fname2,fnameout):
    embed_ratios(fname1)
    embed_ratios(fname2)

    f17 = r.TFile(fname1)
    f18 = r.TFile(fname2)

    b_eff_17 = f17.Get("h2_BTaggingEff_csv_med_Eff_b")
    b_eff_18 = f18.Get("h2_BTaggingEff_csv_med_Eff_b")
    c_eff_17 = f17.Get("h2_BTaggingEff_csv_med_Eff_c")
    c_eff_18 = f18.Get("h2_BTaggingEff_csv_med_Eff_c")
    l_eff_17 = f17.Get("h2_BTaggingEff_csv_med_Eff_udsg")
    l_eff_18 = f18.Get("h2_BTaggingEff_csv_med_Eff_udsg")

    fratio = r.TFile(fnameout, "recreate")

    b_eff_ratio = b_eff_18.Clone("b_eff_ratio_18over17")
    b_eff_ratio.Divide(b_eff_17)

    c_eff_ratio = c_eff_18.Clone("c_eff_ratio_18over17")
    c_eff_ratio.Divide(c_eff_17)

    l_eff_ratio = l_eff_18.Clone("l_eff_ratio_18over17")
    l_eff_ratio.Divide(l_eff_17)

    b_eff_ratio.Write()
    c_eff_ratio.Write()
    l_eff_ratio.Write()

    fratio.Write()
    fratio.Close()

if __name__ == "__main__":

    # embed_ratios("btageff__ttbar_amc_102X_deepCSV.root")
    # embed_ratios("btageff__SMS-T1tttt_2016_80X_deepCSV.root")
    # embed_ratios("btageff__SMS-T1tttt_2017_94X_deepCSV.root")
    # embed_ratios("btageff__qcd_2017_94X_deepCSV.root")

    # make_two_year_ratio(
    #         "output_tt_2017.root",
    #         "output_tt_2018.root",
    #         "hists_bcl_eff_ratios_18over17.root",
    #         )

    make_two_year_ratio(
            "btageff__qcd_2017_94X_deepCSV.root",
            "btageff__qcd_2018_102X_deepCSV.root",
            "hists_qcd_bcl_eff_ratios_18over17.root",
            )


