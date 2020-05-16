import os
import sys
import glob
import analysis.utils.pyrun as pyrun
import ROOT as r
import numpy as np
r.gROOT.SetBatch()

def make_chunks(x, size=5):
    return map(list,np.array_split(x,len(x)//size))

def run_chain((index,info)):
    ch, options, outname = info
    ret = r.ScanChain(ch, options, outname)
    return index, ret

if __name__ == "__main__":
    r.gROOT.ProcessLine(".L ../CORE/CMS3_CORE.so")
    r.gROOT.ProcessLine(".L ScanChain.C+")

    ch = r.TChain("Events")

    # options = "2018,quiet"
    # fnames = glob.glob("/hadoop/cms/store/group/snt/run2_mc2018//TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2_MINIAODSIM_CMS4_V10-02-04/merged_ntuple_1*.root")
    # outputdir = "outputs_2018"
    # outputname = "btageff__ttbar_amc_102X_deepCSV.root"

    # options = "2017,quiet"
    # fnames = glob.glob("/hadoop/cms/store/group/snt/run2_mc2017//TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05/merged_ntuple_1*.root")
    # outputdir = "outputs_2017"
    # outputname = "btageff__ttbar_amc_94X_deepCSV.root"

    # options = "2016,80X,quiet"
    # fnames = glob.glob("/hadoop/cms/store/user/namin/run2_moriond17_cms4/ProjectMetis/SMS-T1tttt_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2-PUSpring16Fast_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1_MINIAODSIM_CMS4_V00-00-02_2017Sep27/merged_ntuple*.root")
    # outputdir = "outputs_2016"
    # outputname = "btageff__SMS-T1tttt_2016_80X_deepCSV.root"

    # options = "2017,quiet"
    # fnames = glob.glob("/hadoop/cms/store/group/snt/run2_mc2017//SMS-T1tttt_TuneCP2_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PUFall17Fast_pilot_94X_mc2017_realistic_v15-v3_MINIAODSIM_CMS4_V10-02-05/merged_ntuple_1*.root")
    # outputdir = "outputs_2017"
    # outputname = "btageff__SMS-T1tttt_2017_94X_deepCSV.root"

    # options = "2017,quiet"
    # fnames = glob.glob("/hadoop/cms/store/group/snt/run2_mc2017//QCD_Pt_80to120_TuneCP5_13TeV_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_MINIAODSIM_CMS4_V10-02-05//*.root")
    # outputdir = "outputs_2017_qcd"
    # outputname = "btageff__qcd_2017_94X_deepCSV.root"

    options = "2018,quiet"
    fnames = glob.glob("/hadoop/cms/store/group/snt/run2_mc2018//QCD_Pt_80to120_TuneCP5_13TeV_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04/*.root")
    outputdir = "outputs_2018_qcd"
    outputname = "btageff__qcd_2018_102X_deepCSV.root"

    chunks = make_chunks(fnames,size=1)

    os.nice(10)

    to_run = []
    for ichunk,chunk in enumerate(chunks):
        ch = r.TChain("Events")
        for fname in chunk:
            ch.Add(fname)
        outname = "{}/output_{}_of_{}.root".format(outputdir,ichunk+1,len(chunks))
        to_run.append([ch,options,outname])

    os.system("mkdir -p {}".format(outputdir))

    run_chain(next(enumerate(to_run)))

    # print to_run
    # map(run_chain, enumerate(to_run))

    # Delete older files
    os.system("rm {}/*.root".format(outputdir))

    runner = pyrun.Runner(nproc=min(len(to_run),15), func=run_chain, dot_type=2)
    runner.add_args(to_run)
    runner.run()

    # Hadd to get final file
    os.system("hadd -k -f {} {}/*.root".format(outputname,outputdir))
