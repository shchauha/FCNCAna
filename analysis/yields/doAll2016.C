{
    gROOT->ProcessLine(".L ../misc/class_files/v8.02/SS.cc+");
    gROOT->ProcessLine(".L ../../common/CORE/Tools/dorky/dorky.cc+");
    gROOT->ProcessLine(".L ScanChain.C+");
    
    //TString basedir_data = "/home/users/namin/2018/fourtop/all/FTAnalysis/analysis/yields/local/year_2016/";
    //TString basedir_data = "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24/output/year_2016/";
    TString basedir_data = "/home/shchauha/2019/Analysis/Babies/year_2016/skimfix2/";
    //TString basedir_data = "/home/shchauha/2019/Analysis/Babies/year_2016/";
    TString basedir_mc = basedir_data;

    /* Options:
     *     useInclusiveSFs: use inclusive (ie not era-dependent) lepton scale factors
     *     doFakes: Run data-driven fake estimate
     *     doFlips: Run charge mis-id estimate
     *     zeroMissingInnerHits: Require exactly zero missing inner hits
     *     evaluateBDT: calculates event BDT, other event_bdt==0
     */
    TString options = "useInclusiveSFs Data2016 doSS";
    //TString outputdir = "v4.2_data";
    TString outputdir = "v4.2_skimfix2";
    //TString outputdir = "v";
    gSystem->Exec(Form("rm -rf %s", outputdir.Data()));
    gSystem->Exec(Form("mkdir -p %s", outputdir.Data()));

    // gSystem->Exec(Form("git show HEAD > %s/git.diff", outputdir.Data()));
    // gSystem->Exec(Form("git diff >> %s/git.diff", outputdir.Data()));

    //Data
    TChain ch_data("t", "data");
    ch_data.Add(basedir_data+"Data*.root");
    ScanChain(&ch_data, options, outputdir);

    // // Data-Driven Fakes
    // TChain ch_fakes("t", "fakes");
    // ch_fakes.Add(basedir_data+"Data*.root");
    // ScanChain(&ch_fakes, options + "doFakes", outputdir);

    // // Data-Driven Flips
    // TChain ch_flips("t", "flips");
    // ch_flips.Add(basedir_data+"Data*.root");
    // ScanChain(&ch_flips, options + "doFlips", outputdir);

   
    // TChain ch_fcnc("t", "fcnc");
    // ch_fcnc.Add("TT_FCNC_hut_top.root");
    // ScanChain(&ch_fcnc, options, outputdir);
    

    // // Monte-Carlo Backgrounds
    // TChain ch_ttw("t", "ttw");
    // ch_ttw.Add(basedir_mc+"TTWnlo.root");
    // ScanChain(&ch_ttw, options, outputdir);

    // TChain ch_ttz("t", "ttz");
    // ch_ttz.Add(basedir_mc+"TTZnlo.root");
    // ScanChain(&ch_ttz, options, outputdir);

    // TChain ch_tth("t", "tth");
    // ch_tth.Add(basedir_mc+"TTHtoNonBB.root");
    // ScanChain(&ch_tth, options, outputdir);
    
    // TChain ch_dy("t", "dy");
    // //ch_dy.Add(basedir_mc+"DY_low.root");
    // ch_dy.Add(basedir_mc+"DY_high.root");
    // ScanChain(&ch_dy, options, outputdir);

    // TChain ch_wjets("t", "wjets");
    // ch_wjets.Add(basedir_mc+"WJets*.root");
    // ScanChain(&ch_wjets, options, outputdir);

    // // TChain ch_tt("t", "tt");
    // // ch_tt.Add(basedir_mc+"TTBAR*.root");
    // // ScanChain(&ch_tt, options, outputdir);

    // // TChain ch_vv("t", "vv");
    // // ch_vv.Add(basedir_mc+"WZ.root");
    // // /* ch_vv.Add(basedir_mc+"WW.root"); */
    // // ch_vv.Add(basedir_mc+"WWDPS.root");
    // // ch_vv.Add(basedir_mc+"ZZ.root");
    // // ScanChain(&ch_vv, options, outputdir);

    // TChain ch_wz("t", "wz");
    // ch_wz.Add(basedir_mc+"WZ.root");
    // ScanChain(&ch_wz, options, outputdir);

    // TChain ch_ww("t", "ww");
    // ch_ww.Add(basedir_mc+"WW.root");
    // ScanChain(&ch_ww, options, outputdir);

    // TChain ch_rares("t", "rares");
    // ch_rares.Add(basedir_mc+"WZG.root");
    // ch_rares.Add(basedir_mc+"WZZ.root");
    // ch_rares.Add(basedir_mc+"ZZZ.root");
    // ch_rares.Add(basedir_mc+"WWZ.root");
    // ch_rares.Add(basedir_mc+"WWW.root");
    // ch_rares.Add(basedir_mc+"TZQ.root");
    // ch_rares.Add(basedir_mc+"TTZZ.root");
    // ch_rares.Add(basedir_mc+"TTWZ.root");
    // ch_rares.Add(basedir_mc+"TTWW.root");
    // ch_rares.Add(basedir_mc+"TTWH.root");
    // ch_rares.Add(basedir_mc+"TTHH.root");
    // ch_rares.Add(basedir_mc+"TTTW.root");
    // ch_rares.Add(basedir_mc+"TTTJ.root");
    // ScanChain(&ch_rares, options, outputdir);

    // TChain ch_singletop("t", "singletop");
    // ch_singletop.Add(basedir_mc+"ST1.root");
    // ch_singletop.Add(basedir_mc+"ST2.root");
    // ScanChain(&ch_singletop, options, outputdir);

    // TChain ch_ttdl("t", "ttdl");
    // ch_ttdl.Add(basedir_mc+"TTDL.root");
    // ScanChain(&ch_ttdl, options, outputdir);

    // TChain ch_fcnc("t", "fcnc");
    // ch_fcnc.Add(basedir_mc+"FCNC_hut_top.root");
    // ch_fcnc.Add(basedir_mc+"FCNC_hut_antitop.root");
    // ScanChain(&ch_fcnc, options, outputdir);

    // TChain ch_ttsl("t", "ttsl");
    // ch_ttsl.Add(basedir_mc+"TTSLtop.root");
    // ch_ttsl.Add(basedir_mc+"TTSLtopbar.root");
    // ScanChain(&ch_ttsl, options, outputdir);

}

