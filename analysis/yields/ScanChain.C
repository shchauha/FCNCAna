#pragma GCC diagnostic ignored "-Wsign-compare"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"
#include "TMVA/Reader.h"

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

#include "../misc/class_files/v8.02/SS.h"
#include "../../common/CORE/Tools/dorky/dorky.h"
#include "../../common/CORE/Tools/utils.h"
#include "../misc/common_utils.h"
#include "../misc/bdt.h"

#include "../misc/tqdm.h"

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> Vec4;
// using namespace tas;

// float lumiAG = getLumi();
bool STOP_REQUESTED = false;
// float lumiAG = 36.3;
struct HistCol2D {
    map<string, TH2D> in;
    HistCol2D(vector<string> regions, const string& name, int nbinsx, float lowx, float highx, int nbinsy, float lowy, float highy, vector<HistCol2D*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            in.emplace(region, TH2D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbinsx, lowx, highx, nbinsy, lowy, highy));
        }
        if (registry != nullptr) registry->push_back(this);
    }
    void Fill(const string& region, int id1, int id2, float valx, float valy, float weight) { in[region].Fill(valx, valy, weight); }
    void Write() { for (auto p : in) p.second.Write(); }
};
struct HistCol {
  map<string, TH1D> in;


    HistCol(vector<string> regions, const string& name, int nbins, const float* bins, vector<HistCol*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, bins));
        }
        if (registry != nullptr)
            registry->push_back(this);
    }

    HistCol(vector<string> regions, const string& name, int nbins, float low, float high, vector<HistCol*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, low, high));
        }
        if (registry != nullptr)
            registry->push_back(this);
    }

    void Fill(const string& region, int id1, int id2, float val, float weight) {
        in[region].Fill(val, weight);
    }

    void Write() {
      for (auto p : in) p.second.Write();
    }
};


float calcDeltaPhi(float phi1, float phi2){
  float dPhi = phi1 - phi2;
  while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
  return fabs(dPhi);
}

float calcMT(float pt1, float phi1, float pt2, float phi2){
  return sqrt( 2 * pt1 * pt2 * ( 1 - cos( phi1 - phi2 ) ) );
}


float isr_reweight(bool useIsrWeight, int year, int nisrmatch) {
    if (!useIsrWeight) return 1;
    if (ss::is_real_data()) return 1;
    return isrWeight(year, nisrmatch, 10); // 10 is ttbar
}

float nb_reweight(int nbtags) {
    if (ss::is_real_data()) return 1;
    std::vector<float> weights = { 1.06, 0.96, 0.99, 1.29, 1.31 };
    // std::vector<float> weights = { 1.00, 1.00, 1.00, 1.29, 1.00 }; // FIXME only reweighting nb==3
    return weights[min(nbtags,4)];
}


int ScanChain(TChain *ch, TString options="", TString outputdir="outputs"){
  
  signal(SIGINT, [](int){
      cout << "SIGINT Caught, stopping after current event" << endl;
      STOP_REQUESTED=true;
    });
    STOP_REQUESTED=false;
    
    bool doFakes = options.Contains("doFakes");
    bool doTTHF = options.Contains("doTTHF");
    bool doNonTTHF = options.Contains("doNonTTHF");
    bool doFlips = options.Contains("doFlips");
    bool useInclusiveSFs = options.Contains("useInclusiveSFs");
    bool zeroMissingInnerHits = options.Contains("zeroMissingInnerHits");
    bool evaluateBDT = options.Contains("evaluateBDT");
    bool useIsoTriggers = options.Contains("useIsoTriggers");
    bool useNonIsoTriggers = options.Contains("useNonIsoTriggers");
    bool doHighHT = options.Contains("doHighHT");
    bool doSS = options.Contains("doSS");
    bool doHEMBefore = options.Contains("doHEMBefore");
    bool doHEMAfter = options.Contains("doHEMAfter");
    bool noLeptonPtCut = options.Contains("noLeptonPtCut");
    bool doTruthFake = options.Contains("doTruthFake");
    bool useNewMET = options.Contains("useNewMET");
    bool quiet = options.Contains("quiet");
    bool minPtFake18 = options.Contains("minPtFake18");
    bool new2016FRBins = options.Contains("new2016FRBins");
    bool noISRWeights = options.Contains("noISRWeights");
    ana_t analysis = FTANA;
    if (doSS) {
        analysis = SSANA;
    }

    TString proc(ch->GetTitle());
    // bool useIsrWeight = proc.Contains("tt_");
    bool useIsrWeight = proc.Contains("tt_") 
        or proc.Contains("ttw_") 
        or proc.Contains("ttz_") 
        or proc.Contains("tth_")
        or proc.Contains("tthf_")
        or proc.Contains("ttnonhf_");

    bool doScaleUnc = proc.Contains("tt_");
    bool useTTBB = proc.Contains("tt_") 
        or proc.Contains("ttw_") 
        or proc.Contains("ttz_") 
        or proc.Contains("tth_");

    if (noISRWeights) useIsrWeight = false;

    // We may have derived the fake rate map throwing away leptons with pT<18 (e.g., 2017), so
    // we need to apply this cut here to be consistent
    // float min_pt_fake = minPtFake18 ? 18. : -1;
    float min_pt_fake = minPtFake18 ? 18. : -1;

    int year;
    float lumiAG = 0.;
    bool is2016(false), is2017(false), is2018(false);
    if (options.Contains("Data2016")) {
        lumiAG = getLumi(2016);
        year = 2016;
        is2016 = true;
    } else if (options.Contains("Data2017")) {
        lumiAG = useNonIsoTriggers ? 36.529: getLumi(2017);
        year = 2017;
        is2017 = true;
    } else if (options.Contains("Data2018")) {
        lumiAG = getLumi(2018);
        year = 2018;
        is2018 = true;
    } else {
        cout << "Which Need to specify year!\n";
        return -1;
    }
    // Clear already-seen list
    duplicate_removal::clear_list();

    // Used to determine which "era" a MC event is in
    TRandom *tr1 = new TRandom();

    if (!quiet) cout << "Working on " << proc << endl;


    vector<string> regions = {

      "br",                          // OS tight-tight and variants
      "sshh",                        // HH SS
      
    };
    // doHighHT = true; // FIXME FIXME

    vector<HistCol*> registry;
    vector<HistCol2D*> registry2D;
    // HistCol h_met         (regions, "met"        , 30, 0   , 300 , &registry);
    HistCol h_met         (regions, "met"        , 60, 0   , 600 , &registry);
    HistCol h_metphi      (regions, "metphi"     , 60, -3.2, 3.2 , &registry);
    HistCol h_rawmet      (regions, "rawmet"     , 60, 0   , 150 , &registry);
    //HistCol h_calomet     (regions, "calomet"    , 60, 0   , 600 , &registry);
    HistCol h_ht          (regions, "ht"         , 16, 0   , 1600, &registry);
    HistCol h_htb         (regions, "htb"        , 16, 0   , 1600, &registry);
    HistCol h_mll         (regions, "mll"        , 30, 0   , 300 , &registry);
    HistCol h_mllzoom     (regions, "mllzoom"    , 60, 60  , 120 , &registry);
    HistCol h_zmll        (regions, "zmll"       , 30, 0   , 300 , &registry);
    HistCol h_m3l        (regions, "m3l"       , 20, 60   , 120 , &registry);
    HistCol h_mtmin       (regions, "mtmin"      , 15, 0   , 300 , &registry);
    HistCol h_dphil1met      (regions, "dphil1met"     , 60, -3.2, 3.2 , &registry);
    HistCol h_dphil2met      (regions, "dphil2met"     , 60, -3.2, 3.2 , &registry);
    HistCol h_dphimetj1      (regions, "dphimetj1"     , 60, -3.2, 3.2 , &registry);
    HistCol h_metmmht      (regions, "metmmht"     , 60, 0, 2.0 , &registry);
    HistCol h_metoverptj1      (regions, "metoverptj1"     , 50, 0.,2., &registry);
    HistCol h_nleps       (regions, "nleps"      , 5, -0.5 , 4.5 , &registry);
    HistCol h_njets       (regions, "njets"      , 8 , 0   , 8   , &registry);
    HistCol h_nisrjets    (regions, "nisrjets"   , 5 , 0   , 5   , &registry);
    HistCol h_nisrmatch   (regions, "nisrmatch"  , 5 , 0   , 5   , &registry);
    HistCol h_nlb40       (regions, "nlb40"      , 5 , 0   , 5   , &registry);
    HistCol h_ntb40       (regions, "ntb40"      , 8 , 0   , 8   , &registry);
    HistCol h_nbtags      (regions, "nbtags"     , 5 , 0   , 5   , &registry);
    HistCol h_nbtagsheavyup      (regions, "nbtagsheavyup"     , 5 , 0   , 5   , &registry);
    HistCol h_nbtagsheavydown      (regions, "nbtagsheavydown"     , 5 , 0   , 5   , &registry);
    HistCol h_nbtagslightup      (regions, "nbtagslightup"     , 5 , 0   , 5   , &registry);
    HistCol h_nbtagslightdown      (regions, "nbtagslightdown"     , 5 , 0   , 5   , &registry);
    HistCol h_nbnj      (regions, "nbnj"     , 20 , 0   , 20   , &registry);
    HistCol h_bdisc1      (regions, "bdisc1"     , 100,0.4 , 1.0 , &registry);
    HistCol h_maxmjoverpt (regions, "maxmjoverpt", 50, 0   , 0.35, &registry);
    HistCol h_btagid      (regions, "btagid"     , 100 , 0   , 100   , &registry);

    HistCol h_pt1         (regions, "pt1"        , 30, 0   , 300 , &registry);
    HistCol h_pt2         (regions, "pt2"        , 30, 0   , 300 , &registry);
    HistCol h_pt3         (regions, "pt3"        , 30, 0   , 300 , &registry);
    HistCol h_pte         (regions, "pte"        , 30, 0   , 300 , &registry);
    HistCol h_ptm         (regions, "ptm"        , 30, 0   , 300 , &registry);

    HistCol h_eta1        (regions, "eta1"       , 25, -3.2, 3.2 , &registry);
    HistCol h_eta2        (regions, "eta2"       , 25, -3.2, 3.2 , &registry);
    HistCol h_etae        (regions, "etae"       , 25, -3.2, 3.2 , &registry);
    HistCol h_etam        (regions, "etam"       , 25, -3.2, 3.2 , &registry);
    HistCol h_etaelnt        (regions, "etaelnt"       , 25, -3.2, 3.2 , &registry);
    HistCol h_etamlnt        (regions, "etamlnt"       , 25, -3.2, 3.2 , &registry);

    HistCol h_phie        (regions, "phie"       , 50, -3.2, 3.2 , &registry);
    HistCol h_phim        (regions, "phim"       , 50, -3.2, 3.2 , &registry);
    HistCol h_cphie        (regions, "cphie"       , 15, -3.2, 3.2 , &registry);
    HistCol h_cphim        (regions, "cphim"       , 15, -3.2, 3.2 , &registry);

    HistCol h_ptrel1      (regions, "ptrel1"     , 15, 0   , 50  , &registry);
    HistCol h_ptrel2      (regions, "ptrel2"     , 15, 0   , 50  , &registry);
    HistCol h_ptrellnt    (regions, "ptrellnt"   , 15, 0   , 50  , &registry);
    HistCol h_ptrele      (regions, "ptrele"     , 15, 0   , 50  , &registry);
    HistCol h_ptrelm      (regions, "ptrelm"     , 15, 0   , 50  , &registry);

    HistCol h_ptrelfail1  (regions, "ptrelfail1" , 15, 0   , 50  , &registry);
    HistCol h_ptrelfail2  (regions, "ptrelfail2" , 15, 0   , 50  , &registry);
    HistCol h_ptrelfaillnt(regions, "ptrelfaillnt", 15, 0   , 50  , &registry);
    HistCol h_ptrelfaile  (regions, "ptrelfaile" , 15, 0   , 50  , &registry);
    HistCol h_ptrelfailm  (regions, "ptrelfailm" , 15, 0   , 50  , &registry);

    HistCol h_ptratio1    (regions, "ptratio1"   , 30, 0   , 1.5 , &registry);
    HistCol h_ptratio2    (regions, "ptratio2"   , 30, 0   , 1.5 , &registry);
    HistCol h_ptratiolnt  (regions, "ptratiolnt" , 30, 0   , 1.5 , &registry);
    HistCol h_ptratioe    (regions, "ptratioe"   , 30, 0   , 1.5 , &registry);
    HistCol h_ptratiom    (regions, "ptratiom"   , 30, 0   , 1.5 , &registry);

    HistCol h_miniiso1    (regions, "miniiso1"   , 15,  0  , 0.2 , &registry);
    HistCol h_miniiso2    (regions, "miniiso2"   , 15,  0  , 0.2 , &registry);
    HistCol h_miniisolnt  (regions, "miniisolnt" , 15,  0  , 0.2 , &registry);
    HistCol h_miniisoe    (regions, "miniisoe"   , 15,  0  , 0.2 , &registry);
    HistCol h_miniisom    (regions, "miniisom"   , 15,  0  , 0.2 , &registry);

    HistCol h_nmiss1      (regions, "nmiss1"     , 3 , -0.5, 2.5 , &registry);
    HistCol h_nmiss2      (regions, "nmiss2"     , 3 , -0.5, 2.5 , &registry);

    HistCol h_dphil1l2    (regions, "dphil1l2"   , 15,  0  , 4   , &registry);
    HistCol h_detal1l2    (regions, "detal1l2"   , 30,  -4 , 4   , &registry);
    HistCol h_absdetal1l2 (regions, "absdetal1l2", 15,  0  , 4   , &registry);

    HistCol h_type        (regions, "type"       , 3 , -0.5, 2.5 , &registry);
    HistCol h_q1          (regions, "q1"         , 2 , -2  , 2   , &registry);
    HistCol h_type3l      (regions, "type3l"     , 4 , -0.5, 3.5 , &registry);
    HistCol h_nvtx        (regions, "nvtx"       , 31, -0.5, 61.5, &registry);

    // HistCol h_ptj1        (regions, "ptj1"       , 50, 0   , 500 , &registry);
    HistCol h_ptj1        (regions, "ptj1"       , 50, 0   , 1000 , &registry);
    HistCol h_ptj2        (regions, "ptj2"       , 50, 0   , 500 , &registry);
    HistCol h_ptj3        (regions, "ptj3"       , 50, 0   , 500 , &registry);
    HistCol h_ptj4        (regions, "ptj4"       , 50, 0   , 500 , &registry);
    HistCol h_ptj5        (regions, "ptj5"       , 50, 0   , 500 , &registry);
    HistCol h_ptj6        (regions, "ptj6"       , 30, 0   , 300 , &registry);
    HistCol h_ptj7        (regions, "ptj7"       , 30, 0   , 300 , &registry);
    HistCol h_ptj8        (regions, "ptj8"       , 30, 0   , 300 , &registry);

    HistCol h_avgcdisc (regions          , "avgcdisc"           , 30, -2 , 1 , &registry);
    HistCol h_nforwardjets20 (regions    , "nforwardjets20"     , 10, 0 , 10 , &registry);
    HistCol h_ntrijets (regions          , "ntrijets"           , 30, 0 , 30 , &registry);
    HistCol h_trijet_meandisc (regions   , "trijet_meandisc"    , 30, 0 , 1 , &registry);
    HistCol h_trijet_leadingdisc (regions, "trijet_leadingdisc" , 30, 0 , 1 , &registry);
    HistCol h_trijet_subleadingdisc (regions, "trijet_subleadingdisc" , 30, 0 , 1 , &registry);
    HistCol h_trijet_numhigh (regions    , "trijet_numhigh"     , 10, 0 , 10 , &registry);
    HistCol h_trijet_frachigh (regions   , "trijet_frachigh"    , 30, 0 , 1 , &registry);

    HistCol h_ml1j1       (regions, "ml1j1"      , 30, 0   , 300 , &registry);
    HistCol h_matchtype   (regions, "matchtype"  , 4 , -0.5, 3.5 , &registry);

    HistCol2D h_ptabsetae       (regions, "ptetae"      , 40,0,400,30,0.,3., &registry2D);
    HistCol2D h_ptabsetam       (regions, "ptetam"      , 40,0,400,30,0.,3., &registry2D);


    // Declare a bunch of event variables to be filled below in the loop
    float lep1ccpt, lep2ccpt, lep3ccpt;
    float lep1pt,   lep2pt,   lep3pt;
    float lep1eta,  lep2eta,  lep3eta;
    float lep1phi,  lep2phi,  lep3phi;
    int   lep1id,   lep2id,   lep3id;
    int   lep1q,    lep2q,    lep3q;
    int   lep1good, lep2good, lep3good;

    float lep1ptrel, lep2ptrel;
    float lep1miniiso, lep2miniiso;
    float lep1ptratio, lep2ptratio;
    float dphil1l2, detal1l2;
    float nmiss1, nmiss2;

    int nleps, njets, nbtags;
    float ht, htb, met, metphi, rawmet, calomet;
    float bdisc1;
    float maxmjoverpt, ml1j1;
    int matchtype;
    int nlb40, ntb40, nisrjets, nisrmatch;
    float ptj1, ptj2, ptj3, ptj4, ptj5, ptj6, ptj7, ptj8;
    float jet1phi;

    int SR;
    float weight;    
    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();

    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);

    tqdm bar;
    bar.set_theme_braille();

    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        if (STOP_REQUESTED) break;

        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("t");
        samesign.Init(tree);

        TString filename(currentFile->GetTitle());

        auto tokens = filename.Tokenize("/");
        auto basename = ((TObjString*)(tokens->At(tokens->GetEntries()-1)))->String().Data();
        bar.set_label(basename);

        for(unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

            if (STOP_REQUESTED) break;

            samesign.GetEntry(event);
            nEventsTotal++;
            
            if (!quiet) bar.progress(nEventsTotal, nEventsChain);

            // Simple cuts first to speed things up
            lep1ccpt = ss::lep1_coneCorrPt();
            lep2ccpt = ss::lep2_coneCorrPt();

            if (not noLeptonPtCut) {
                if (lep1ccpt < 25) continue;
                if (lep2ccpt < 20) continue;
            }
            lep1id = ss::lep1_id();
            lep2id = ss::lep2_id();


            bool pass_trig = (doSS) ? ss::fired_trigger_ss() : ss::fired_trigger();
	    
            if (!pass_trig) continue;
            if (!ss::passes_met_filters()) continue;
            // Reject duplicates
            if (ss::is_real_data()) {
                duplicate_removal::DorkyEventIdentifier id(ss::run(), ss::event(), ss::lumi());
                if (duplicate_removal::is_duplicate(id)) continue;
            }

            // Save a bunch of event info for quick reference later
            njets = ss::njets();
            nbtags = ss::nbtags();
	    bdisc1 = 1;
            met = ss::met();
            metphi = ss::metPhi();
            rawmet = ss::rawmet();
            lep3ccpt = ss::lep3_coneCorrPt();
            lep1pt = ss::lep1_p4().pt();
            lep2pt = ss::lep2_p4().pt();
            lep3pt = ss::lep3_p4().pt();
            lep1eta = ss::lep1_p4().eta();
            lep2eta = ss::lep2_p4().eta();
            lep3eta = ss::lep3_p4().eta();
            lep1phi = ss::lep1_p4().phi();
            lep2phi = ss::lep2_p4().phi();
            lep3phi = ss::lep3_p4().phi();
            lep3id = ss::lep3_id();
            lep1q = ss::bdt_q1();
            lep2q = (lep2id > 0) ? -1 : 1;
            lep3q = (lep3id > 0) ? -1 : 1;
            lep1good = ss::lep1_passes_id();
            lep2good = ss::lep2_passes_id();
            lep3good = ss::lep3_passes_id();
            lep1ptrel = ss::lep1_ptrel_v1();
            lep2ptrel = ss::lep2_ptrel_v1();
            lep1miniiso = ss::lep1_miniIso();
            lep2miniiso = ss::lep2_miniIso();
            nmiss1 = ss::lep1_el_exp_innerlayers();
            nmiss2 = ss::lep2_el_exp_innerlayers();	   
	    lep1ptratio = ss::lep1_ptratio();
	    lep2ptratio = ss::lep2_ptratio();
            
            if (doSS) {
                nleps = (lep3good) ? ((ss::lep4_passes_id() and (ss::lep4_p4().pt() > (abs(ss::lep4_id())==11 ? 15 : 10))) ? 4 : 3) : 2;
            } else {
                nleps = (lep3good and lep3ccpt > 20) ? ((ss::lep4_passes_id() and ss::lep4_p4().pt() > 20) ? 4 : 3) : 2;
            }
            ht = ss::ht();
            nisrmatch = ss::nisrMatch();
            SR = 0; // Just a dummy for now
	    
            /* hyp_class
             * 1: SS, loose-loose
             * 2: SS, tight-loose (or loose-tight)
             * 3: SS, tight-tight
             * 4: OS, tight-tight
             * 5: SS, inSituFR
             * 6: SS, tight-tight and fails Z-veto (lies! hyp_class==6 != lep1good and lep2good)
             */
            int hyp_class = ss::hyp_class();



            //Calculate weight
            weight = ss::is_real_data() ? 1 : ss::scale1fb()*lumiAG;

            if (!ss::is_real_data()) {
                weight *= getTruePUw(year, ss::trueNumInt()[0]);
                if (lep1good) weight *= leptonScaleFactor(year, lep1id, lep1ccpt, lep1eta, ht);
                if (lep2good) weight *= leptonScaleFactor(year, lep2id, lep2ccpt, lep2eta, ht);
                if (not doSS) {
                    if (lep3good) weight *= leptonScaleFactor(year, lep3id, lep3ccpt, lep3eta, ht);
                }
                if (doSS or !lep3good) {
                    weight *= triggerScaleFactor(year, lep1id, lep2id, lep1pt, lep2pt, lep1eta, lep2eta, ht, analysis, 0);
                }
                weight *= ss::weight_btagsf();
                weight *= isr_reweight(useIsrWeight, year, nisrmatch);
                if (year == 2016) weight *= ss::prefire2016_sf();
                if (year == 2017) weight *= ss::prefire2017_sf();
		
                if (useTTBB and (ss::extragenb() >= 2)) weight *= 1.7;
                
            }


            bool class6Fake = false;
            if (doFakes) {
                if (hyp_class == 6) {
                    bool lep1_lowpt_veto = lep1pt < (abs(lep1id) == 11 ? 15 : 10);
                    bool lep2_lowpt_veto = lep2pt < (abs(lep2id) == 11 ? 15 : 10);
                    bool lep3_lowpt_veto = lep3pt < (abs(lep3id) == 11 ? 15 : 10);
                    int nfakes = 0;
                    if (ss::lep3_fo() and !ss::lep3_tight() and !lep3_lowpt_veto and lep1good and lep2good && lep3pt>min_pt_fake) {  // lep3 fake
                        float fr = fakeRate(year, lep3id, lep3ccpt, lep3eta, ht, analysis, new2016FRBins, !minPtFake18);
                        class6Fake = true;
                        nfakes++;
                        weight *= fr / (1-fr);
                    }
                    if (ss::lep2_fo() and !ss::lep2_tight() and !lep2_lowpt_veto and lep1good and lep3good && lep2pt>min_pt_fake) {  // lep2 fake
                        float fr = fakeRate(year, lep2id, lep2ccpt, lep2eta, ht, analysis, new2016FRBins, !minPtFake18);
                        class6Fake = true;
                        nfakes++;
                        weight *= fr / (1-fr);
                    }
                    if (ss::lep1_fo() and !ss::lep1_tight() and !lep1_lowpt_veto and lep2good and lep3good && lep1pt>min_pt_fake) {  // lep1 fake
                        float fr = fakeRate(year, lep1id, lep1ccpt, lep1eta, ht, analysis, new2016FRBins, !minPtFake18);
                        class6Fake = true;
                        nfakes++;
                        weight *= fr / (1-fr);
                    }
                    if (!class6Fake) {
                        continue; // No fakes!
                    }
                    if (nfakes == 2) weight *= -1;
                } else if (hyp_class == 1 or hyp_class == 2) {
                    bool foundGoodLoose = false;
                    if (ss::lep1_passes_id()==0 && lep1pt>min_pt_fake) {
                        float fr = fakeRate(year, lep1id, lep1ccpt, lep1eta, ht, analysis, new2016FRBins, !minPtFake18);
                        weight *= fr/(1.-fr);
                        foundGoodLoose = true;
                    }
                    if (ss::lep2_passes_id()==0 && lep2pt>min_pt_fake) {
                        float fr = fakeRate(year, lep2id, lep2ccpt, lep2eta, ht, analysis, new2016FRBins, !minPtFake18);
                        weight *= fr/(1.-fr);
                        foundGoodLoose = true;
                    }
                    if (!foundGoodLoose)
                        continue;
                    // subtract double FO (why is this?)
                    if (hyp_class == 1 && lep1pt>min_pt_fake && lep2pt>min_pt_fake) weight *= -1.;
                    hyp_class = 3; // we've faked a SS Tight-Tight with a SS LL or SS TL
                                   // Basically just update this so it gets put in the SR
                } else {
                    continue; // Not a fakeing hyp_class
                }
            }


            if (doFlips) {
                if (hyp_class == 4) hyp_class = 3; // we've flipped an OS to a SS
                // else if (hyp_class == 6) class6Fake = true;
                else continue;
                float flipFact = 0.;
                if (abs(lep1id) == 11) {
                    float flr = flipRate(year, lep1pt, lep1eta);
                    flipFact += (flr/(1-flr));
                }
                if (abs(lep2id) == 11) {
                    float flr = flipRate(year, lep2pt, lep2eta);
                    flipFact += (flr/(1-flr));
                }
                weight *= flipFact;
                if (weight == 0.0) continue; // just quit if there are no flips.
            }

            if (not doSS) {
                // if all 3 charges are the same, throw the event away
                if (nleps > 2 and ((lep1id>0 and lep2id>0 and lep3id>0) or
                            (lep1id<0 and lep2id<0 and lep3id<0))) continue;
            }

            auto getmll = [](const Vec4& p1, const Vec4& p2, float ccpt1=-1, float ccpt2=-1) {
	      /* Calculate dilepton mass with optional rescaling based on cone-corrected lepton pt */
	      if (ccpt1 == -1) return (p1 + p2).M();
	      else             return (p1*ccpt1/p1.pt() + p2*ccpt2/p2.pt()).M();
            };
            float m12 = getmll(ss::lep1_p4(), ss::lep2_p4(), lep1ccpt, lep2ccpt);
            float m13 = getmll(ss::lep1_p4(), ss::lep3_p4(), lep1ccpt, lep3ccpt);
            float m23 = getmll(ss::lep2_p4(), ss::lep3_p4(), lep2ccpt, lep3ccpt);
            float m3l = (ss::lep1_p4()+ss::lep2_p4()+ss::lep3_p4()).M();

            auto z_cand = [](int id1, int id2, float mll) {
                return abs(id1) == abs(id2) and  // Same flavor
                       id1*id2<0 and             // Opposite sign
                       abs(mll - 91.2) < 15;     // Z-mass window
            };
            bool zcand12 = z_cand(lep1id, lep2id, m12);
            bool zcand13 = z_cand(lep1id, lep3id, m13);
            bool zcand23 = z_cand(lep2id, lep3id, m23);
            float mllos = fabs(m13 - 91.2) < fabs(m23 - 91.2) ? m13 : m23;
	    


	    
            auto fill_region = [&](const string& region, float weight) {
	      if (std::find(regions.begin(), regions.end(), region) == regions.end()) return;
	      
                // Fill all observables for a region
                auto do_fill = [region, lep1id, lep2id, weight](HistCol& h, float val, float extraweight=1.) {
                    h.Fill(region, lep1id, lep2id, val, weight*extraweight);
                };
                auto do_fill2D = [region, lep1id, lep2id, weight](HistCol2D& h, float valx, float valy) {
                    h.Fill(region, lep1id, lep2id, valx, valy, weight);
                };
	        do_fill(h_met, met);
                do_fill(h_metphi, metphi);
                do_fill(h_rawmet, rawmet);
                //do_fill(h_calomet, calomet);
                do_fill(h_ht, ht);
                do_fill(h_htb, htb);
                do_fill(h_mll, m12);
                do_fill(h_mllzoom, m12);
                if (nleps > 2) do_fill(h_zmll, mllos);
                if (nleps > 2) do_fill(h_m3l, m3l);
                do_fill(h_nleps, nleps);
                do_fill(h_mtmin, ss::mtmin());
                do_fill(h_dphil1met, calcDeltaPhi(lep1phi,metphi));
                do_fill(h_dphil2met, calcDeltaPhi(lep2phi,metphi));
                do_fill(h_njets, njets);
                do_fill(h_nisrjets, nisrjets);
                do_fill(h_nisrmatch, nisrmatch);	        
                do_fill(h_ntb40, ntb40);
                do_fill(h_nbtags, nbtags);
                int nbnj = 5*min(nbtags,3)+(max(min(njets,6),2)-2);
                do_fill(h_nbnj, nbnj);

                do_fill(h_pt1, lep1ccpt);
                do_fill(h_pt2, lep2ccpt);
	        if (nleps > 2) do_fill(h_pt3, lep3pt);
                do_fill(h_eta1,   ss::lep1_p4().eta());
	        do_fill(h_eta2,   ss::lep2_p4().eta());

                int looseleg = -1;
                if (hyp_class == 2) {
                    looseleg = (lep1good ? 2 : 1);
                }

                if (looseleg == 1) do_fill(abs(lep1id) == 11 ? h_etaelnt     : h_etamlnt,     lep1eta);
                if (looseleg == 2) do_fill(abs(lep2id) == 11 ? h_etaelnt     : h_etamlnt,     lep2eta);

                do_fill(h_ptrel1, lep1ptrel);
                do_fill(h_ptrel2, lep2ptrel);
                if (looseleg > 0) do_fill(h_ptrellnt, looseleg == 1 ? lep1ptrel   : lep2ptrel);
                do_fill(abs(lep1id) == 11 ? h_ptrele   : h_ptrelm,   lep1ptrel);
                do_fill(abs(lep2id) == 11 ? h_ptrele   : h_ptrelm,   lep2ptrel);

                do_fill2D(abs(lep1id) == 11 ? h_ptabsetae     : h_ptabsetam,     lep1pt, fabs(lep1eta));
                do_fill2D(abs(lep2id) == 11 ? h_ptabsetae     : h_ptabsetam,     lep2pt, fabs(lep2eta));


                do_fill(h_ptratio1, lep1ptratio);
                do_fill(h_ptratio2, lep2ptratio);
	        int type = ss::hyp_type();
                do_fill(h_type,   type>1 ? type-1 : type);
                do_fill(h_q1,   lep1id>0 ? -1 : 1);
                if (nleps > 2) {
                    if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 39) do_fill(h_type3l, 0); // mu mu mu
                    if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 37) do_fill(h_type3l, 1); // mu mu e
                    if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 35) do_fill(h_type3l, 2); // mu e e
                    if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 33) do_fill(h_type3l, 3); // e e e
                }
                do_fill(h_nvtx,   ss::nGoodVertices());
	    
	    };

            bool BR_lite = ht > 300       and njets >= 2 and
                           nbtags >= 2    and met >= 50;
            bool BR = BR_lite and hyp_class == 3;
            if (BR) fill_region("br", weight);            
            if (hyp_class == 3 and njets >= 2 and met > 50
                    and lep1ccpt > 25 and lep2ccpt > 25
                    ) {
                fill_region("sshh", weight);
            }


        }//event loop
        delete file;
    }//file loop

    TFile f1(Form("%s/histos_%s.root", outputdir.Data(), ch->GetTitle()), "RECREATE");
    for (HistCol* coll : registry) coll->Write();
    for (HistCol2D* coll : registry2D) coll->Write();
    f1.Close();
    if (!quiet) cout << "\n Done!" << endl;
    return 0;
}

