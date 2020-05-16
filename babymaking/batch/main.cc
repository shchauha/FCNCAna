#include "./CORE/CMS3.h"
#include "./CORE/SSSelections.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TBranchElement.h"
#include "TString.h"
#include "helper_babymaker.h"
#include "CORE/Tools/goodrun.h"
#include "CORE/Tools/JetCorrector.h"
#include "CORE/Tools/datasetinfo/getDatasetInfo.h"
#include "CORE/Tools/jetcorr/JetCorrectionUncertainty.h"  
#include "CORE/Tools/jetcorr/SimpleJetCorrectionUncertainty.h"
#include "cxxopts.h"

bool STOP_REQUESTED;

int condor_chirp(std::string key, std::string value);

int main(int argc, char *argv[]){


    bool verbose = false;
    bool arg_notfastsim = false;
    bool arg_ignorescale1fb = false;
    bool arg_skipos = false;
    bool arg_xrootd = false;
    bool arg_ignorebadfiles = false;
    auto arg_files = std::vector<std::string>();
    std::string arg_output = "output.root";
    int arg_nevents = -1;
    int arg_year = -1;
    int arg_eventsel = 0;
    std::vector<TString> vfilenames;


    try {

        cxxopts::Options options(argv[0], "FTAnalysis babymaker");
        options.positional_help("[optional args]").show_positional_help();
        options.add_options()
            ("h,help", "Print help")
            ("v,verbose", "Enable verbosity", cxxopts::value<bool>(verbose))
            ("notfastsim", "Don't treat this as fastsim even if it is", cxxopts::value<bool>(arg_notfastsim))
            ("ignorescale1fb", "Ignore the scale1fb.txt file", cxxopts::value<bool>(arg_ignorescale1fb))
            ("ignorebadfiles", "Ignore/skip files with I/O errors", cxxopts::value<bool>(arg_ignorebadfiles))
            ("skipos", "Skip OS events for MC", cxxopts::value<bool>(arg_skipos))
            ("xrootd", "Use xrootd by force", cxxopts::value<bool>(arg_xrootd))
            ("o,output", "Output name", cxxopts::value<std::string>(arg_output)->default_value(arg_output))
            ("n,nevents", "Number of events", cxxopts::value<int>(arg_nevents)->default_value(std::to_string(arg_nevents)))
            ("y,year", "Processing year", cxxopts::value<int>(arg_year)->default_value(std::to_string(arg_year)))
            ("eventsel", "Select particular event by number", cxxopts::value<int>(arg_eventsel)->default_value(std::to_string(arg_eventsel)))
            ("f,files", "Files or file patterns", cxxopts::value<std::vector<std::string>>(arg_files))
            ;
        options.parse_positional({"files"});
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            return 0;
        }

        if (result.count("files")) {
            auto& ff = result["files"].as<std::vector<std::string>>();
            for (const auto& f : ff) {
                vfilenames.push_back(TString(f));
            }
        } else {
            std::cout << ">>> Need a list of files!" << std::endl;
            std::cout << options.help({""}) << std::endl;
            return 1;
        }

        // std::cout <<  " vfilenames.size(): " << vfilenames.size() <<  std::endl;
        // std::cout <<  " verbose: " << verbose <<  std::endl;
        // std::cout <<  " arg_notfastsim: " << arg_notfastsim <<  std::endl;
        // std::cout <<  " arg_ignorescale1fb: " << arg_ignorescale1fb <<  std::endl;
        // std::cout <<  " arg_skipos: " << arg_skipos <<  std::endl;
        // std::cout <<  " arg_output: " << arg_output <<  std::endl;
        // std::cout <<  " arg_nevents: " << arg_nevents <<  std::endl;

    } catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        return 1;
    }


    // std::cout << result << std::endl;

    signal(SIGINT, [](int){
            cout << "SIGINT Caught, stopping after current event" << endl;
            STOP_REQUESTED=true;
            });
    STOP_REQUESTED = false;

    TString filenames = vfilenames.at(0);
    TString outname(arg_output);
    unsigned int nevents_max = (arg_nevents > 0 ? arg_nevents : 0);
    // bool useXrootd = false;
    char hostnamestupid[100];
    int res = gethostname(hostnamestupid, 100);
    TString hostname(hostnamestupid);
    std::cout << ">>> Hostname is " << hostname << std::endl;  
    // bool useXrootd = !(hostname.Contains("t2.ucsd.edu")) or arg_xrootd;
    bool useXrootd = true; // NOTE this is only ~5-10% slower locally, but not a burden on fuse mount!!
    bool goodRunsOnly = true;
    // bool useXrootd = false;
    if (hostname.Contains("uafino")) {
        std::cout << ">>> We're on uafino, so using xrootd!" << std::endl;  
        useXrootd = true;
    }

    std::cout << ">>> Args: " << std::endl;
    std::cout << "     filenames:   " << std::endl;
    for (auto f : vfilenames) { 
        std::cout << "         " << f << std::endl;
    }
    std::cout << "     outname:     " << outname << std::endl;
    std::cout << "     nevents_max: " << nevents_max << std::endl;

    bool isData = filenames.Contains("run2_data") || filenames.Contains("Run201");
    int iSignal = 0; 
    bool isFastsim = 0; 
    int year = -1;
    if (
            filenames.Contains("Run2016")
            || filenames.Contains("Moriond17")
            || filenames.Contains("RunIISummer16")
            || filenames.Contains("run2_data2016")
            || filenames.Contains("run2_moriond17")
            || filenames.Contains("2016_mc_rpv")
            || filenames.Contains("run2_mc2016")
       ) year = 2016;
    if (
            filenames.Contains("Run2017")
            || filenames.Contains("RunIIFall17")
            || filenames.Contains("_mc2017_")
            || filenames.Contains("run2_mc2017")
       ) year = 2017;
    if (
            filenames.Contains("Run2018")
            || filenames.Contains("RunIISpring18")
            || filenames.Contains("RunIISummer18")
            || filenames.Contains("Autumn18")
            || filenames.Contains("run2_mc2018")
       ) year = 2018;

    if (arg_year >= 2016) {
        year = arg_year;
        std::cout << ">>> [!] You manually overrode the year to " << year << std::endl;
    }

    if (year < 0) {
        std::cout << ">>> [!] Couldn't figure out year, so setting it to 2017. Make sure this is what you want!" << std::endl;
        year = 2017;
    } else {
        std::cout << ">>> Figured out that the year is " << year << "." << std::endl;
    }

    // 2016

    string good_run_file;
    string jecEra; // XXX NOTE that you also have to change the run-dependent data JECs below too
    std::string jecEraMC;

    if (year == 2016) {
        gconf.year = year;
        gconf.cmssw_ver = 80;
        gconf.ea_version = 1;
        gconf.btag_disc_wp = 0.6324;
        gconf.WP_DEEPCSV_TIGHT  = 0.8958;
        gconf.WP_DEEPCSV_MEDIUM = 0.6324;
        gconf.WP_DEEPCSV_LOOSE  = 0.2219;
        gconf.multiiso_el_minireliso = 0.12;
        gconf.multiiso_el_ptratio = 0.80;
        gconf.multiiso_el_ptrel = 7.2;
        gconf.multiiso_mu_minireliso = 0.16;
        gconf.multiiso_mu_ptratio = 0.76;
        gconf.multiiso_mu_ptrel = 7.2;
        jecEra = "Summer16_23Sep2016BCDV4";
        jecEraMC = "Summer16_23Sep2016V4";
        if (
                filenames.Contains("17Jul2018") ||
                filenames.Contains("RunIISummer16MiniAODv3-PUMoriond17_94X") ||
                filenames.Contains("94X") ||
                filenames.Contains("MiniAODv3") ||
                filenames.Contains("run2_mc2016_94x")
                ) {
            gconf.cmssw_ver = 94;
            // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
            // Uh, yeah... They changed in the 4th decimal only
            gconf.WP_DEEPCSV_TIGHT  = 0.8953;
            gconf.WP_DEEPCSV_MEDIUM = 0.6321;
            gconf.WP_DEEPCSV_LOOSE  = 0.2217;
            jecEra = "Summer16_07Aug2017BCD_V11";
            jecEraMC = "Summer16_07Aug2017_V11";
        }
        good_run_file = "goodRunList/goldenJson_2016rereco_36p46ifb.txt";
    }
    if (year == 2017) {
        gconf.year = year;
        gconf.ea_version = 4;
        gconf.cmssw_ver = 94;

        gconf.btag_disc_wp = 0.4941;
        gconf.WP_DEEPCSV_TIGHT  = 0.8001;
        gconf.WP_DEEPCSV_MEDIUM = 0.4941;
        gconf.WP_DEEPCSV_LOOSE  = 0.1522;

        // gconf.btag_disc_wp = 0.3033;
        // gconf.WP_DEEPCSV_TIGHT  = 0.7489;
        // gconf.WP_DEEPCSV_MEDIUM = 0.3033;
        // gconf.WP_DEEPCSV_LOOSE  = 0.0521;

        // gconf.multiiso_el_minireliso = 0.09;
        // gconf.multiiso_el_ptratio = 0.85;
        // gconf.multiiso_el_ptrel = 9.2;
        // gconf.multiiso_mu_minireliso = 0.12;
        // gconf.multiiso_mu_ptratio = 0.80;
        // gconf.multiiso_mu_ptrel = 7.5;

        // From FT_jecprefire_18Dec6.pdf
        gconf.multiiso_el_minireliso = 0.07;
        gconf.multiiso_el_ptratio = 0.78;
        gconf.multiiso_el_ptrel = 8.0;
        gconf.multiiso_mu_minireliso = 0.11;
        gconf.multiiso_mu_ptratio = 0.74;
        gconf.multiiso_mu_ptrel = 6.8;

        good_run_file = "goodRunList/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1_snt.txt";
        jecEra = "Fall17_17Nov2017B_V32";
        jecEraMC = "Fall17_17Nov2017_V32";
        // jecEra = "Fall17_17Nov2017B_V6";
        // jecEraMC = "Fall17_17Nov2017_V6";
        gconf.SS_innerlayers = 0;
    }
    if (year == 2018) {
        gconf.year = year;
        gconf.ea_version = 4;
        gconf.cmssw_ver = 102;

        gconf.btag_disc_wp = 0.4184;
        gconf.WP_DEEPCSV_TIGHT  = 0.7527;
        gconf.WP_DEEPCSV_MEDIUM = 0.4184;
        gconf.WP_DEEPCSV_LOOSE  = 0.1241;

        // From FT_jecprefire_18Dec6.pdf
        gconf.multiiso_el_minireliso = 0.07;
        gconf.multiiso_el_ptratio = 0.78;
        gconf.multiiso_el_ptrel = 8.0;
        gconf.multiiso_mu_minireliso = 0.11;
        gconf.multiiso_mu_ptratio = 0.74;
        gconf.multiiso_mu_ptrel = 6.8;

        // good_run_file = "goodRunList/Cert_314472-324209_13TeV_PromptReco_Collisions18_JSON_snt.txt"; // 50.98
        // good_run_file = "goodRunList/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON_snt.txt"; // 58.83
        good_run_file = "goodRunList/goldenJson_2018_final_59p76ifb_snt.txt"; // 59.76ifb (rerecoABCD+promptD + recovered 50ifb + latest preapproved lumi normtag)
        // jecEra = "Fall17_17Nov2017C_V32";
        // jecEraMC = "Fall17_17Nov2017_V32";
        // jecEra = "Autumn18_V3_OldRes";
        // jecEraMC = "Autumn18_V3";
        jecEra = "Autumn18_RunA_V8";
        // jecEra = "Autumn18_RunA_V15"; // Test to compare against V8. Results still use V8
        jecEraMC = "Autumn18_V8";
        // jecEra = "Autumn18_RunA_V5";
        // jecEraMC = "Autumn18_V5";
        gconf.SS_innerlayers = 0;
    }

    if (filenames.Contains("/SMS")) {
        if (filenames.Contains("SMS-T1tttt_")) iSignal = 1;
        else if (filenames.Contains("SMS-T5qqqqVV_")) {
            iSignal = 2;
            if (filenames.Contains("_dM20_")) iSignal = 3;
        }
        else if (filenames.Contains("SMS-T6ttWW")) iSignal = 10;
        else if (filenames.Contains("SMS-T6ttHZ")) iSignal = 11;
        else if (filenames.Contains("SMS-T1ttbb_Tune")) iSignal = 6;
        else if (filenames.Contains("SMS-T1qqqqL")) iSignal = 7;
        else if (filenames.Contains("SMS-T5tttt_dM175_Tune")) iSignal = 5;
        else if (filenames.Contains("SMS-T5ttcc")) iSignal = 4;
        else if (filenames.Contains("SMS-T1tbs")) iSignal = -1; // set later
        else {
            std::cout << ">>> [!] Can't figure out which signal sample this is!" << std::endl;
            return 1;
        }
        isFastsim = true;
    }
    // if (filenames.Contains("_HToTT_")) iSignal = 101; // iSignal > 100 used only for non SMS stuff
    if (filenames.Contains("_HToTT_") && filenames.Contains("Summer16Fast")) {
        iSignal = 101; // fastsim higgs vs fullsim higgs
        isFastsim = true;
    }
    if (filenames.Contains("T1qqqqL") and not filenames.Contains("Fast")) {
        iSignal = 201; // RPV sample
        isFastsim = false;
    }
    if (filenames.Contains("T1tbs")) {
        iSignal = 301; // RPV sample
        isFastsim = false;
    }

    if (arg_notfastsim) {
        iSignal = 0;
        isFastsim = false;
        std::cout << ">>> [!] Not actually a signal file!" << std::endl;
    }

    //Set up file and tree
    babyMaker *mylooper = new babyMaker();
    mylooper->MakeBabyNtuple(outname.Data(), isFastsim, iSignal);


    if (verbose) {
        mylooper->verbose = true;
        std::cout << ">>> [!] Verbose mode" << std::endl;
    }
    if (arg_ignorescale1fb) {
        mylooper->ignore_scale1fb = true;
        std::cout << ">>> [!] Ignoring scale1fb text file!" << std::endl;
    }
    if (arg_skipos) {
        mylooper->ignore_os = true;
        std::cout << ">>> [!] Skipping OS events for MC!" << std::endl;
    }
    if (arg_eventsel > 0) {
        mylooper->evt_cut = arg_eventsel;
        std::cout << ">>> [!] Only processing single event with evt_event==" << arg_eventsel << std::endl;
    }

    // 2018 sample has PSweights but doesn't have the string in the dataset name :(
    if (filenames.Contains("/TTTT") && (filenames.Contains("PSweights") || year == 2018)) {
        mylooper->read_psweights = true;
        std::cout << ">>> [!] This is a TTTT sample and PSweights exist, so reading them" << std::endl;
    }



    TChain *chain = new TChain("Events");

    if (useXrootd) {
        std::cout << ">>> Using xrootd" << std::endl;  
    }

    for (auto fname : vfilenames) {
        if (useXrootd) {
            // fname.ReplaceAll("/hadoop/cms", "root://cmsxrootd.fnal.gov/");
            fname.ReplaceAll("/hadoop/cms", "root://redirector.t2.ucsd.edu/");
        }
        auto nfiles = chain->Add(fname);
        std::cout << ">>> Adding " << nfiles << " file(s): " << fname << std::endl;  
    }

    TObjArray *listOfBranches = chain->GetListOfBranches();
    TIter branchIter(listOfBranches);
    TBranchElement *currentBranch = 0;
    while ( (currentBranch = (TBranchElement*)branchIter.Next()) ) { 
        TString bname = TString(currentBranch->GetName());
        if (bname.Contains("genMaker_genHEPMCweight") and iSignal==0) {
            mylooper->has_lhe_branches = true;
            std::cout << ">>> [!] This sample has the new gen_LHE_* branches, so using them" << std::endl;
            break;
        }
    }

    if (mylooper->has_lhe_branches) {
        if (filenames.Contains("LOoblique")) {
            mylooper->has_lhe_branches = false;
            std::cout << ">>> [!] This sample is private and new gen_LHE_* branches are 0, so falling back to the old method" << std::endl;
        }

    }

    //Add good run list
    set_goodrun_file(good_run_file.c_str());


    //Make Jet Energy Uncertainties
    JetCorrectionUncertainty *jecUnc = 0;
    if (!isData)      jecUnc = new JetCorrectionUncertainty("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_Uncertainty_AK4PFchs.txt"); 
    if (iSignal > 0 and year == 2016) jecUnc = new JetCorrectionUncertainty("CORE/Tools/jetcorr/data/run2_25ns/Spring16_FastSimV1/Spring16_FastSimV1_Uncertainty_AK4PFchs.txt"); 
    if (iSignal > 0 and year == 2017) jecUnc = new JetCorrectionUncertainty("CORE/Tools/jetcorr/data/run2_25ns/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC_Uncertainty_AK4PFchs.txt"); 
    if (iSignal > 0 and year == 2018) jecUnc = new JetCorrectionUncertainty("CORE/Tools/jetcorr/data/run2_25ns/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC_Uncertainty_AK4PFchs.txt"); 

    //Init MVA
    createAndInitMVA("./CORE", true, true, 80);

    //JEC files -- 25 ns MC
    std::vector<std::string> jetcorr_filenames_25ns_MC_pfL1;
    std::vector<std::string> jetcorr_filenames_25ns_MC_pfL2L3;
    std::vector<std::string> jetcorr_filenames_25ns_MC_pfL1L2L3;
    jetcorr_filenames_25ns_MC_pfL1.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L1FastJet_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L2Relative_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L3Absolute_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L1FastJet_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L2Relative_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L3Absolute_AK4PFchs.txt");

    //JEC files -- 25 ns DATA
    std::vector<std::string> jetcorr_filenames_25ns_DATA_pfL1;
    std::vector<std::string> jetcorr_filenames_25ns_DATA_pfL2L3;
    std::vector<std::string> jetcorr_filenames_25ns_DATA_pfL1L2L3;

    //JEC files -- 25 ns FASTSIM
    std::vector<std::string> jetcorr_filenames_25ns_FASTSIM_pfL1;
    std::vector<std::string> jetcorr_filenames_25ns_FASTSIM_pfL2L3;
    std::vector<std::string> jetcorr_filenames_25ns_FASTSIM_pfL1L2L3;
    if (year == 2016) {
        jetcorr_filenames_25ns_FASTSIM_pfL1.push_back("CORE/Tools/jetcorr/data/run2_25ns/Spring16_FastSimV1/Spring16_FastSimV1_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_25ns_FASTSIM_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Spring16_FastSimV1/Spring16_FastSimV1_L2Relative_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Spring16_FastSimV1/Spring16_FastSimV1_L3Absolute_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Spring16_FastSimV1/Spring16_FastSimV1_L2Relative_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Spring16_FastSimV1/Spring16_FastSimV1_L3Absolute_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Spring16_FastSimV1/Spring16_FastSimV1_L1FastJet_AK4PFchs.txt");
    } else if (year == 2017) {
        jetcorr_filenames_25ns_FASTSIM_pfL1.push_back("CORE/Tools/jetcorr/data/run2_25ns/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_25ns_FASTSIM_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC_L2Relative_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC_L3Absolute_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC_L2Relative_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC_L3Absolute_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC_L1FastJet_AK4PFchs.txt");
    } else if (year == 2018) {
        jetcorr_filenames_25ns_FASTSIM_pfL1.push_back("CORE/Tools/jetcorr/data/run2_25ns/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC_L1FastJet_AK4PFchs.txt");
        jetcorr_filenames_25ns_FASTSIM_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC_L2Relative_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC_L3Absolute_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC_L2Relative_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC_L3Absolute_AK4PFchs.txt"); 
        jetcorr_filenames_25ns_FASTSIM_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC_L1FastJet_AK4PFchs.txt");
    }

    //JECs

    FactorizedJetCorrector *jetCorrAG = nullptr;
    FactorizedJetCorrector *jetCorrAG_L1 = nullptr;
    FactorizedJetCorrector *jetCorrAG_L2L3 = nullptr;
    if (!isData) {
        if (!iSignal) {
            jetCorrAG_L1 = makeJetCorrector(jetcorr_filenames_25ns_MC_pfL1);
            jetCorrAG_L2L3 = makeJetCorrector(jetcorr_filenames_25ns_MC_pfL2L3);
            jetCorrAG = makeJetCorrector(jetcorr_filenames_25ns_MC_pfL1L2L3);
        } else {
            jetCorrAG_L1 = makeJetCorrector(jetcorr_filenames_25ns_FASTSIM_pfL1);
            jetCorrAG_L2L3 = makeJetCorrector(jetcorr_filenames_25ns_FASTSIM_pfL2L3);
            jetCorrAG = makeJetCorrector(jetcorr_filenames_25ns_FASTSIM_pfL1L2L3);
        }
        // for closestJet() in CORE/IsolationTools.cc
        gconf.jet_corrector_L1 = jetCorrAG_L1;
        gconf.jet_corrector_L2L3 = jetCorrAG_L2L3;
    }

    TH2D* count_hist = new TH2D("counts","",560,0,2800,500,0,2500);
    TH1D* weight_hist = new TH1D("weight","",1200,0,1200);

    auto t0 = std::chrono::steady_clock::now();
    auto tlast = std::chrono::steady_clock::now();

    // Set up loop over chain
    unsigned int nEvents = chain->GetEntries();
    unsigned int nEventsTotal = 0;
    unsigned int nEventsToDo = chain->GetEntries();
    if(nEvents > 0) nEventsToDo = nEvents;
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;

    // Event Counting
    std::cout << ">>> Events in TChain: " << nEvents << std::endl;  


    // File Loop
    while ( (currentFile = (TFile*)fileIter.Next()) ) { 
        if (STOP_REQUESTED) break;

        // Get File Content
        if(nEventsTotal >= nEventsToDo) continue;
        if ((nevents_max > 0) && (nEventsTotal >= nevents_max)) break;
        TFile *file = TFile::Open( currentFile->GetTitle() );
        if (!file and arg_ignorebadfiles) {
            std::cout << "Ignore bad file option is set, so continuing onto the next file" << std::endl;
            continue;
        }
        TTree *tree = (TTree*)file->Get("Events");
        TString filename(currentFile->GetTitle());
        cms3.Init(tree);

        // Loop over Events in current file
        unsigned int nEventsTree = tree->GetEntriesFast();
        for (unsigned int evt = 0; evt < nEventsTree; evt++){

            if (STOP_REQUESTED) break;

            if ((nevents_max > 0) && (evt > nevents_max)) break;

            //Get Event Content
            cms3.GetEntry(evt);
            chain->LoadTree(evt);

            nEventsTotal++;

            // HAVE TO MAKE DATA JEC HERE SINCE WE NEED RUN NUMBER BECAUSE SO MANY JECS
            // Only need to check first event because a merged file can't span eras
            if (isData && evt == 0) {
                if (     tas::evt_run() <= 276811) jecEra = (gconf.cmssw_ver == 80 ? "Summer16_23Sep2016BCDV4" : "Summer16_07Aug2017BCD_V11");
                else if (tas::evt_run() <= 278801 && tas::evt_run() >= 276831) jecEra = (gconf.cmssw_ver == 80 ? "Summer16_23Sep2016EFV4" : "Summer16_07Aug2017EF_V11");
                else if (tas::evt_run() <= 280385 && tas::evt_run() >= 278802) jecEra = (gconf.cmssw_ver == 80 ? "Summer16_23Sep2016GV4" : "Summer16_07Aug2017GH_V11");
                else if (tas::evt_run() <  294645 && tas::evt_run() >= 280919) jecEra = (gconf.cmssw_ver == 80 ? "Summer16_23Sep2016HV4" : "Summer16_07Aug2017GH_V11");

                else if (tas::evt_run() <= 299329 && tas::evt_run() >= 297046) jecEra = "Fall17_17Nov2017B_V32";
                else if (tas::evt_run() <= 302029 && tas::evt_run() >= 299368) jecEra = "Fall17_17Nov2017C_V32";
                else if (tas::evt_run() <= 304797 && tas::evt_run() >= 302030) jecEra = "Fall17_17Nov2017DE_V32";
                else if (tas::evt_run() <= 306462 && tas::evt_run() >= 305040) jecEra = "Fall17_17Nov2017F_V32";

                else if (tas::evt_run() <= 316995 && tas::evt_run() >= 315252) jecEra = "Autumn18_RunA_V8";
                else if (tas::evt_run() <= 319312 && tas::evt_run() >= 316998) jecEra = "Autumn18_RunB_V8";
                else if (tas::evt_run() <= 320393 && tas::evt_run() >= 319313) jecEra = "Autumn18_RunC_V8";
                else if (tas::evt_run() <= 325273 && tas::evt_run() >= 320394) jecEra = "Autumn18_RunD_V8";
                // else if (tas::evt_run() <= 316995 && tas::evt_run() >= 315252) jecEra = "Autumn18_RunA_V15"; // See note earlier about V15 vs V8
                // else if (tas::evt_run() <= 319312 && tas::evt_run() >= 316998) jecEra = "Autumn18_RunB_V15";
                // else if (tas::evt_run() <= 320393 && tas::evt_run() >= 319313) jecEra = "Autumn18_RunC_V15";
                // else if (tas::evt_run() <= 325273 && tas::evt_run() >= 320394) jecEra = "Autumn18_RunD_V15";

                else {
                    std::cout << ">>> [!] Shouldn't get here! Can't figure out JEC. isData,run = " << isData << "," << tas::evt_run() << std::endl;
                    break;
                }

                // CORE/Tools/jetcorr/data/run2_25ns/Autumn18_RunB_V5_DATA/Autumn18_RunB_V5_DATA_L2L3Residual_AK4PFchs.txt

                jecUnc = new JetCorrectionUncertainty("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_Uncertainty_AK4PFchs.txt"); 
                jetcorr_filenames_25ns_DATA_pfL1.clear();
                jetcorr_filenames_25ns_DATA_pfL1.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L1FastJet_AK4PFchs.txt");
                jetcorr_filenames_25ns_DATA_pfL2L3.clear();
                jetcorr_filenames_25ns_DATA_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L2Relative_AK4PFchs.txt");
                jetcorr_filenames_25ns_DATA_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L3Absolute_AK4PFchs.txt");
                jetcorr_filenames_25ns_DATA_pfL2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L2L3Residual_AK4PFchs.txt");
                jetcorr_filenames_25ns_DATA_pfL1L2L3.clear();
                jetcorr_filenames_25ns_DATA_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L1FastJet_AK4PFchs.txt");
                jetcorr_filenames_25ns_DATA_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L2Relative_AK4PFchs.txt");
                jetcorr_filenames_25ns_DATA_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L3Absolute_AK4PFchs.txt");
                jetcorr_filenames_25ns_DATA_pfL1L2L3.push_back("CORE/Tools/jetcorr/data/run2_25ns/"+jecEra+"_DATA/"+jecEra+"_DATA_L2L3Residual_AK4PFchs.txt");
                jetCorrAG_L1 = makeJetCorrector(jetcorr_filenames_25ns_DATA_pfL1);
                jetCorrAG_L2L3 = makeJetCorrector(jetcorr_filenames_25ns_DATA_pfL2L3);
                jetCorrAG = makeJetCorrector(jetcorr_filenames_25ns_DATA_pfL1L2L3);
                // for closestJet() in CORE/IsolationTools.cc
                gconf.jet_corrector_L1 = jetCorrAG_L1;
                gconf.jet_corrector_L2L3 = jetCorrAG_L2L3;

            }

            //If data, check good run list
            if (goodRunsOnly && isData && tas::evt_isRealData() && !goodrun(tas::evt_run(), tas::evt_lumiBlock())) continue;

            if (mylooper->evt_cut>0 && tas::evt_event() != mylooper->evt_cut) continue;

            //Progress bar
            CMS3::progress(nEventsTotal, nEvents);

            if (nEventsTotal % 5000 == 0) {
                auto t1 = std::chrono::steady_clock::now();
                float duration = 0.001*(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - tlast)).count();
                float total_duration = 0.001*(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)).count();
                std::cout << "Begin processing entry " << nEventsTotal << " at " << time(0) << " (" << 0.001*round(5000/duration) << " kHz)." << std::endl;
                tlast = t1;
                if (total_duration > 1.5*3600 and nEventsTotal % 50000 == 0) {
                    condor_chirp("ChirpCMSSWEvents",std::to_string((int)nEventsTotal));
                    condor_chirp("ChirpCMSSWElapsed",std::to_string((int)total_duration));
                    condor_chirp("ChirpCMSSWLastUpdate",std::to_string((int)(time(0))));
                    condor_chirp("ChirpMetisStatus",std::string("\"running\""));
                    condor_chirp("ChirpMetisExpectedNevents",std::to_string((int)nEventsToDo));
                    condor_chirp("ChirpCMSSWEventRate",std::to_string(round(5000/duration))); // keep EventRate at end so we know when chirp train is done
                }
            }

            // csErr_t csErr = mylooper->ProcessBaby(filename.Data(), jetCorrAG, jecUnc, iSignal); 
            try {
                mylooper->ProcessBaby(filename.Data(), jetCorrAG, jecUnc, isFastsim, iSignal); 
            } catch (const std::ios_base::failure& e) {
                std::cout << "Got an I/O error!" << std::endl;
                std::cout << e.what() << std::endl;
                if (arg_ignorebadfiles) {
                    std::cout << "Ignore bad file option is set, so continuing onto the next file" << std::endl;
                    break;
                }
            }
            // int SR = csErr.SR; 
            // bool isGood = csErr.isGood; 

            //c-s error variables
            if (isFastsim) {
                vector<float> sparms = tas::sparm_values();
                float s1 = 1.;
                float s2 = 1.;
                if (sparms.size() > 0) s1 = sparms[0];
                if (sparms.size() > 1) s2 = sparms[1];
                count_hist->Fill(s1,s2);
            } else {
                count_hist->Fill(1,1);
            }

            if (!tas::evt_isRealData()) {

                if (!mylooper->has_lhe_branches) {

                    auto genweights = tas::genweights();
                    if (genweights.size()>110) {
                        float nom = genweights[0];
                        float scale_up_raw = genweights[4];
                        float scale_down_raw = genweights[8];
                        if (isFastsim) {
                            // fastsim shifted by 1
                            nom = genweights[1];
                            scale_up_raw = genweights[5];
                            scale_down_raw = genweights[9];
                        }
                        float sum_pdf = 0.;
                        float sum2_pdf = 0.;
                        int N = 100;
                        for (int ipdf = 9; ipdf < 9+N; ipdf++) {
                            sum_pdf += genweights[ipdf];
                            sum2_pdf += pow(genweights[ipdf],2);
                        }
                        float rms = sqrt(max(sum2_pdf/N - pow(sum_pdf/N,2),(double)0.0));
                        // in 2017, variations are Hessian, so multiply by sqrt(N-1)
                        if (year == 2017) {
                            rms *= sqrt(99);
                        }
                        float pdf_up_raw = (sum_pdf/N+rms);
                        float pdf_down_raw = (sum_pdf/N-rms);
                        if (fabs(pdf_up_raw/nom) > 8.) pdf_up_raw = (pdf_up_raw/nom > 0 ? 8.*nom : -8.*nom);
                        if (fabs(pdf_down_raw/nom) > 8.) pdf_down_raw = (pdf_down_raw/nom > 0 ? 8.*nom : -8.*nom);
                        weight_hist->Fill(0.5, nom);
                        weight_hist->Fill(1.5, scale_up_raw);
                        weight_hist->Fill(2.5, scale_down_raw);
                        weight_hist->Fill(3.5, pdf_up_raw);
                        weight_hist->Fill(4.5, pdf_down_raw);
                    } else {
                        // Fallback, but should be checked always
                        weight_hist->Fill(0.5, 1.);
                        weight_hist->Fill(1.5, 1.);
                        weight_hist->Fill(2.5, 1.);
                        weight_hist->Fill(3.5, 1.);
                        weight_hist->Fill(4.5, 1.);
                    }

                } else {

                        float nom = tas::genHEPMCweight();
                        float scale_up_raw = tas::gen_LHEweight_QCDscale_muR2_muF2()*nom;
                        float scale_down_raw = tas::gen_LHEweight_QCDscale_muR0p5_muF0p5()*nom;
                        float pdf_up_raw = tas::gen_LHEweight_PDFVariation_Up()*nom;
                        float pdf_down_raw = tas::gen_LHEweight_PDFVariation_Dn()*nom;
                        if (fabs(pdf_up_raw/nom) > 8.) pdf_up_raw = (pdf_up_raw/nom > 0 ? 8.*nom : -8.*nom);
                        if (fabs(pdf_down_raw/nom) > 8.) pdf_down_raw = (pdf_down_raw/nom > 0 ? 8.*nom : -8.*nom);
                        weight_hist->Fill(0.5, nom);
                        weight_hist->Fill(1.5, scale_up_raw);
                        weight_hist->Fill(2.5, scale_down_raw);
                        weight_hist->Fill(3.5, pdf_up_raw);
                        weight_hist->Fill(4.5, pdf_down_raw);

                }

            }
            
        } // event loop 

        // // XXX hmm, these 2 lines cause segfaults. why? maybe ROOT garbage collection already
        // // takes care of this when file goes out of scope...
        // file->Close();
        // delete file;

        cout << ">>> File done" <<endl;

    } // file loop

    auto t1 = std::chrono::steady_clock::now();
    float total_duration = 0.001*(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)).count();
    std::cout << "Processed " << nEventsTotal << " events in " << total_duration << " seconds @ " << 0.001*round(nEventsTotal/total_duration) << " kHz" << std::endl;
    std::cout << "Stored " << mylooper->BabyTree->GetEntries() << " events in the output tree." << std::endl;

    if (total_duration > 1.5*3600) {
        condor_chirp("ChirpCMSSWEvents",std::to_string((int)nEventsTotal));
        condor_chirp("ChirpCMSSWElapsed",std::to_string((int)total_duration));
        condor_chirp("ChirpCMSSWLastUpdate",std::to_string((int)(time(0))));
        condor_chirp("ChirpMetisStatus",std::string("\"done\""));
        condor_chirp("ChirpMetisExpectedNevents",std::to_string((int)nEventsToDo));
        condor_chirp("ChirpCMSSWEventRate",std::to_string(round(nEventsTotal/total_duration))); // keep EventRate at end so we know when chirp train is done
    }


    // Cleanup
    mylooper->BabyFile->cd();
    mylooper->BabyTree->Write();
    count_hist->Write();
    weight_hist->Write();
    mylooper->BabyFile->Close();

    if ((nEvents != nEventsTotal) || nEvents==0) {
        std::cout << "Number of input events in tree is " << nEvents << " but we processed " << nEventsTotal << ", ";
        std::cout << "so something went wrong and the return code will not be 0." << std::endl;
        return 1;
    }

    return 0;

}

int condor_chirp(std::string key, std::string value) {
    std::string cmd = "condor_chirp set_attr_delayed " + key + " " + value;
    std::cout << "Running: " << cmd << std::endl;
    system(cmd.c_str());
    return 0;
}

/* vim: set ft=cpp: */
