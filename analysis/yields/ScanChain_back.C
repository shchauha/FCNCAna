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
    map<string, TH1D> ee;
    map<string, TH1D> em;
    map<string, TH1D> mm;
    map<string, TH1D> in;


    HistCol(vector<string> regions, const string& name, int nbins, const float* bins, vector<HistCol*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            ee.emplace(region, TH1D((base_name + "_ee").c_str(), (base_title + " ee").c_str(), nbins, bins));
            em.emplace(region, TH1D((base_name + "_em").c_str(), (base_title + " em").c_str(), nbins, bins));
            mm.emplace(region, TH1D((base_name + "_mm").c_str(), (base_title + " mm").c_str(), nbins, bins));
            in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, bins));
        }
        if (registry != nullptr)
            registry->push_back(this);
    }

    HistCol(vector<string> regions, const string& name, int nbins, float low, float high, vector<HistCol*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            ee.emplace(region, TH1D((base_name + "_ee").c_str(), (base_title + " ee").c_str(), nbins, low, high));
            em.emplace(region, TH1D((base_name + "_em").c_str(), (base_title + " em").c_str(), nbins, low, high));
            mm.emplace(region, TH1D((base_name + "_mm").c_str(), (base_title + " mm").c_str(), nbins, low, high));
            in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, low, high));
        }
        if (registry != nullptr)
            registry->push_back(this);
    }

    void Fill(const string& region, int id1, int id2, float val, float weight) {
        in[region].Fill(val, weight);
        if (abs(id1) == 11 and abs(id2) == 11) {
            ee[region].Fill(val, weight);
        } else if (abs(id1) == 13 and abs(id2) == 13) {
            mm[region].Fill(val, weight);
        } else if ((abs(id1) == 11 and abs(id2) == 13) or
                   (abs(id1) == 13 and abs(id2) == 11)) {
            em[region].Fill(val, weight);
        } else {
            cout << "These ids are garbage: (" << id1 << ", " << id2 << ")\n";
        }
    }

    void Write() {
        for (auto p : ee) p.second.Write();
        for (auto p : em) p.second.Write();
        for (auto p : mm) p.second.Write();
        for (auto p : in) p.second.Write();
    }
};

tuple<int, int, int, float, float> calc_jet_quants(int year) {
    // 2017
    float bloose =  0.1522;
    float bmed = 0.4941;
    float btight = 0.8001;
    if (year == 2016) {
        bloose =  0.2219;
        bmed = 0.6324;
        btight = 0.8958;
    }
    else if (year == 2018) {
        bloose  = 0.1241;
        bmed = 0.4184;
        btight  = 0.7527;
    }
    
    int nisrjets = 0;
    int nlb40 = ss::bdt_nlb40();
    int ntb40 = ss::bdt_ntb40();
    float htb = ss::bdt_htb();
    float maxmjoverpt = ss::bdt_maxmjoverpt();
    for (unsigned int ijet = 0; ijet < ss::jets().size(); ijet++) {
        const auto& jet = ss::jets()[ijet];
        // float disc = ss::jets_disc().at(ijet);
        // if (disc < bmed) nisrjets++;
    }
    return std::make_tuple(nlb40, ntb40, nisrjets, maxmjoverpt, htb);

}

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

float nvtx_reweight_2018(int nvtx) {
    if (nvtx >= -0.5 and nvtx < 1.5) return 15.000;
    if (nvtx >= 1.5 and nvtx < 3.5) return 0.975;
    if (nvtx >= 3.5 and nvtx < 5.5) return 0.893;
    if (nvtx >= 5.5 and nvtx < 7.5) return 0.939;
    if (nvtx >= 7.5 and nvtx < 9.5) return 0.911;
    if (nvtx >= 9.5 and nvtx < 11.5) return 0.886;
    if (nvtx >= 11.5 and nvtx < 13.5) return 0.927;
    if (nvtx >= 13.5 and nvtx < 15.5) return 0.896;
    if (nvtx >= 15.5 and nvtx < 17.5) return 0.929;
    if (nvtx >= 17.5 and nvtx < 19.5) return 0.918;
    if (nvtx >= 19.5 and nvtx < 21.5) return 0.906;
    if (nvtx >= 21.5 and nvtx < 23.5) return 0.912;
    if (nvtx >= 23.5 and nvtx < 25.5) return 0.924;
    if (nvtx >= 25.5 and nvtx < 27.5) return 0.945;
    if (nvtx >= 27.5 and nvtx < 29.5) return 0.946;
    if (nvtx >= 29.5 and nvtx < 31.5) return 1.001;
    if (nvtx >= 31.5 and nvtx < 33.5) return 1.047;
    if (nvtx >= 33.5 and nvtx < 35.5) return 1.130;
    if (nvtx >= 35.5 and nvtx < 37.5) return 1.218;
    if (nvtx >= 37.5 and nvtx < 39.5) return 1.315;
    if (nvtx >= 39.5 and nvtx < 41.5) return 1.518;
    if (nvtx >= 41.5 and nvtx < 43.5) return 1.728;
    if (nvtx >= 43.5 and nvtx < 45.5) return 1.965;
    if (nvtx >= 45.5 and nvtx < 47.5) return 2.325;
    if (nvtx >= 47.5 and nvtx < 49.5) return 2.774;
    if (nvtx >= 49.5 and nvtx < 51.5) return 3.117;
    if (nvtx >= 51.5 and nvtx < 53.5) return 3.877;
    if (nvtx >= 53.5 and nvtx < 55.5) return 4.913;
    if (nvtx >= 55.5 and nvtx < 57.5) return 5.335;
    if (nvtx >= 57.5 and nvtx < 59.5) return 7.496;
    if (nvtx >= 59.5 and nvtx < 66.5) return 14.400;
    return 1.;
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

        "os",                          // OS tight-tight and variants
        // "ftos",                          // OS tight-tight and variants with 25,20 pt
        // "osloose",                    // DY dominated CR
        "htnb1",                    // fake dominated CR
        "tl",                          // SS tight-loose
        // // "osnbrw",                          // OS tight-tight and variants

        // for SS
        // "xgcr",
        "tlnomet",
        "osnomet",
        // "osnomethighmt",
        "lowmetlowht",
        "lowmetallht",
        "tlmet",
        "osmet",
        // "osmethigh",
        // "osmetjet",
        "mlonz",
        // "mlonzhigh",
        // "mlonz0j",

        // "oshtnb3",

        // "tt_isr",                      // TTBar ISR Reweighting derivation region

        // "hhos", "hhosloose", "hhtl",
        // "os_noht",                          // OS tight-tight and variants
        // "br",
        // "crw", "crz",                  // CRZ, CRW
        // "htnb1",                        // fake-dominated CR
        // "htnb1mu",                        // fake-dominated CR
        // "nj2",                        // fake-dominated CR
        // "sshh",                        // HH SS sum

        // "os_btagreweight", // |
        // "os_highbdt", "os_lowbdt",     // |
        // "bdt_nb", "bdt_ht",            // Baseline w/ inverted nbtags/Ht/MET selection
        // "bdt_met","bdt_met_ht",        // |
        // "bdt_train",                   // BDT Training region - BR + CRW
        // "tt_isr_reweight_check",       // |
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

    float bdt_bins[] = {-1.00, -0.48, -0.37, -0.25, -0.14, -0.02, 0.075, 0.21, 0.31, 0.465, 0.60, 0.715, 0.81, 1.00};
    // HistCol h_event_bdt   (regions, "event_bdt"  , 13 , bdt_bins  , &registry);
    // HistCol h_event_bdt   (regions, "eventbdt"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_btagheavyup   (regions, "eventbdt_btagheavyup"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_btagheavydn   (regions, "eventbdt_btagheavydn"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_btaglightup   (regions, "eventbdt_btaglightup"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_btaglightdn   (regions, "eventbdt_btaglightdn"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_scaleup   (regions, "eventbdt_scaleup"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_scaledn   (regions, "eventbdt_scaledn"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_jecup   (regions, "eventbdt_jecup"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_jecdn   (regions, "eventbdt_jecdn"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_jerup   (regions, "eventbdt_jerup"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_jerdn   (regions, "eventbdt_jerdn"  , 15 , 0., 1. , &registry);
    // HistCol h_event_bdt_bbup   (regions, "eventbdt_bbup"  , 15 , 0., 1. , &registry);


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

    // Because TMVA is dumb, we need also float versions of these variables.
    float nbtags_f, njets_f, nlb40_f, ntb40_f, nleps_f, lep1q_f, SR_f;
    auto set_float_vals = [&]() {
        nbtags_f = nbtags;
        njets_f = njets;
        nlb40_f = nlb40;
        ntb40_f = ntb40;
        nleps_f = nleps;
        lep1q_f = lep1q;
        SR_f = SR;
    };
    // TMVA::Reader reader("Silent");
    // reader.AddVariable("nbtags",      &nbtags_f);
    // reader.AddVariable("njets",       &njets_f);
    // reader.AddVariable("met",         &met);
    // reader.AddVariable("ptl2",        &lep2ccpt);
    // reader.AddVariable("nlb40",       &nlb40_f);
    // reader.AddVariable("ntb40",       &ntb40_f);
    // reader.AddVariable("nleps",       &nleps_f);
    // reader.AddVariable("htb",         &htb);
    // reader.AddVariable("q1",          &lep1q_f);
    // reader.AddVariable("ptj1",        &ptj1);
    // reader.AddVariable("ptj6",        &ptj6);
    // reader.AddVariable("ptj7",        &ptj7);
    // reader.AddVariable("ml1j1",       &ml1j1);
    // reader.AddVariable("dphil1l2",    &dphil1l2);
    // reader.AddVariable("maxmjoverpt", &maxmjoverpt);
    // reader.AddVariable("ptl1",        &lep1ccpt);
    // reader.AddVariable("detal1l2",    &detal1l2);
    // reader.AddVariable("ptj8",        &ptj8);
    // reader.AddVariable("ptl3",        &lep3ccpt);

    // reader.AddSpectator("weight",     &weight);
    // reader.AddSpectator("ptl1",       &lep1ccpt);
    // reader.AddSpectator("ptl2",       &lep2ccpt);
    // reader.AddSpectator("SR",         &SR_f);

    //reader.BookMVA("BDT","TMVAClassification_BDT_19vars.xml");


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

            /* if (event > 20) break; */

            // // FIXME
            // if (ss::njets() < 2) continue;
            // if (ss::met() > 50.) continue;
            // if (ss::hyp_class() != 6) continue;

            // Simple cuts first to speed things up
            lep1ccpt = ss::lep1_coneCorrPt();
            lep2ccpt = ss::lep2_coneCorrPt();

            if (not noLeptonPtCut) {
                if (lep1ccpt < 25) continue;
                if (lep2ccpt < 20) continue;
            }

            if (doTTHF and (ss::extragenb() < 1)) continue;
            if (doNonTTHF and (ss::extragenb() >= 1)) continue;

            if (doHEMAfter and ss::run() < 319077) continue; // keep after
            if (doHEMBefore and ss::run() >= 319077) continue; // keep before

            lep1id = ss::lep1_id();
            lep2id = ss::lep2_id();

            if (doHighHT and ss::ht() < 300.) continue;

            // FIXME fix once 2017 ntuples are remade
            // bool pass_trig = (doSS and year != 2017) ? ss::fired_trigger_ss() : ss::fired_trigger();
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
            //bdisc1 = (nbtags >= 1) ? ss::btags_disc()[0] : -1;
	    bdisc1 = 1;
            met = ss::met();
            metphi = ss::metPhi();
            rawmet = ss::rawmet();
            //calomet = ss::calomet();
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
            // if (is2016) {
	    //   lep1ptratio = ss::lep1_closeJet().pt() > 0 ? lep1pt/ss::lep1_closeJet().pt() : 1.0;
	    //   lep2ptratio = ss::lep2_closeJet().pt() > 0 ? lep2pt/ss::lep2_closeJet().pt() : 1.0;
            // }

	   
	    lep1ptratio = ss::lep1_ptratio();
	    lep2ptratio = ss::lep2_ptratio();
            
            if (doSS) {
                nleps = (lep3good) ? ((ss::lep4_passes_id() and (ss::lep4_p4().pt() > (abs(ss::lep4_id())==11 ? 15 : 10))) ? 4 : 3) : 2;
            } else {
                nleps = (lep3good and lep3ccpt > 20) ? ((ss::lep4_passes_id() and ss::lep4_p4().pt() > 20) ? 4 : 3) : 2;
            }
            ht = ss::ht();
            std::tie(nlb40, ntb40, nisrjets, maxmjoverpt, htb) = calc_jet_quants(year);
            nisrmatch = ss::nisrMatch();
            SR = 0; // Just a dummy for now

            // for SS also show 4 leading jets
            // ptj1 = (njets >= 1) ? ss::jets()[0].pt()*ss::jets_undoJEC()[0]*ss::jets_JEC()[0] : 0;
            // ptj2 = (njets >= 2) ? ss::jets()[1].pt()*ss::jets_undoJEC()[1]*ss::jets_JEC()[1] : 0;
            // ptj3 = (njets >= 3) ? ss::jets()[2].pt()*ss::jets_undoJEC()[2]*ss::jets_JEC()[2] : 0;
            // ptj4 = (njets >= 4) ? ss::jets()[3].pt()*ss::jets_undoJEC()[3]*ss::jets_JEC()[3] : 0;
            // jet1phi = (njets >= 1) ? ss::jets()[0].phi() : 0.;
            // if (not doSS) {
            //     ptj5 = (njets >= 5) ? ss::jets()[4].pt()*ss::jets_undoJEC()[4]*ss::jets_JEC()[4] : 0;
            //     ptj6 = (njets >= 6) ? ss::jets()[5].pt()*ss::jets_undoJEC()[5]*ss::jets_JEC()[5] : 0;
            //     ptj7 = (njets >= 7) ? ss::jets()[6].pt()*ss::jets_undoJEC()[6]*ss::jets_JEC()[6] : 0;
            //     ptj8 = (njets >= 8) ? ss::jets()[7].pt()*ss::jets_undoJEC()[7]*ss::jets_JEC()[7] : 0;
            //     ml1j1 = ss::bdt_ml1j1();
            //     set_float_vals();
            // }

            /* hyp_class
             * 1: SS, loose-loose
             * 2: SS, tight-loose (or loose-tight)
             * 3: SS, tight-tight
             * 4: OS, tight-tight
             * 5: SS, inSituFR
             * 6: SS, tight-tight and fails Z-veto (lies! hyp_class==6 != lep1good and lep2good)
             */
            int hyp_class = ss::hyp_class();

            matchtype = 0;
            if (doTruthFake) {
                int nbadlegs = (ss::lep1_motherID() <= 0) + (ss::lep2_motherID() <= 0);
                int ngoodlegs = (ss::lep1_motherID() == 1) + (ss::lep2_motherID() == 1);
                // skip the event if it's truth matched to be prompt prompt.
                // We only want reco tight-tight events that are prompt-nonprompt (or nonprompt nonprompt)
                if (ngoodlegs == 2) continue;

                //consider only prompt or lights
                if ((ss::lep2_motherID()==1 || ss::lep2_motherID()==0) and
                    (ss::lep1_motherID()==1 || ss::lep1_motherID()==0)) matchtype = 0;
                //consider only prompt or bs
                if ((ss::lep2_motherID()==1 || ss::lep2_motherID()==-1) and
                    (ss::lep1_motherID()==1 || ss::lep1_motherID()==-1)) matchtype = 1;
                //consider only prompt or cs
                if ((ss::lep2_motherID()==1 || ss::lep2_motherID()==-2) and
                    (ss::lep1_motherID()==1 || ss::lep1_motherID()==-2)) matchtype = 2;
            }

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

                // // FIXME FIXME this is for testing only
                // if (year == 2018) weight *= nvtx_reweight_2018(ss::nGoodVertices());
            }



            if (doSS and !ss::is_real_data()) {
                if (lep1ccpt < 25. and lep2ccpt < 25. and nleps==2 and year==2017) {
                    // Low low dilepton+HT triggers missing in 2017B (11.6% of 2017 lumi, so weight down MC by (1-0.116));
                    weight *= 1-0.116;
                }
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

            float mvavalue = -1.;
            float mvavalueup = -1.;
            float mvavaluedn = -1.;
            float mvavalueJERup = -1.;
            float mvavalueJERdn = -1.;
            int nrealb = 0;
            if (not doSS) {
                // for (auto flav : ss::btags_flavor()) {
                //     if (abs(flav)==5) nrealb++;
                // }
                float f_nbtags = ss::bdt_nbtags();
                float f_njets = ss::bdt_njets();
                float f_met = ss::bdt_met();
                float f_ptl2 = ss::bdt_ptl2();
                float f_nlb40 = ss::bdt_nlb40();
                float f_ntb40 = ss::bdt_ntb40();
                float f_nleps = ss::bdt_nleps();
                float f_htb = ss::bdt_htb();
                float f_q1 = ss::bdt_q1();
                float f_ptj1 = ss::bdt_ptj1();
                float f_ptj6 = ss::bdt_ptj6();
                float f_ptj7 = ss::bdt_ptj7();
                float f_ml1j1 = ss::bdt_ml1j1();
                float f_dphil1l2 = ss::bdt_dphil1l2();
                float f_maxmjoverpt = ss::bdt_maxmjoverpt();
                float f_ptl1 = ss::bdt_ptl1();
                float f_detal1l2 = ss::bdt_detal1l2();
                float f_ptj8 = ss::bdt_ptj8();
                float f_ptl3 = ss::bdt_ptl3();

                mvavalue = get_prediction(f_nbtags,f_njets,f_met,f_ptl2,f_nlb40,f_ntb40,f_nleps,f_htb,f_q1,f_ptj1,f_ptj6,f_ptj7,f_ml1j1,f_dphil1l2,f_maxmjoverpt,f_ptl1,f_detal1l2,f_ptj8,f_ptl3);
                f_nbtags = ss::bdt_jec_up_nbtags();
                f_njets = ss::bdt_jec_up_njets();
                f_met = ss::bdt_jec_up_met();
                f_htb = ss::bdt_jec_up_htb();
                f_nlb40 = ss::bdt_jec_up_nlb40();
                f_ntb40 = ss::bdt_jec_up_ntb40();
                mvavalueup = get_prediction(f_nbtags,f_njets,f_met,f_ptl2,f_nlb40,f_ntb40,f_nleps,f_htb,f_q1,f_ptj1,f_ptj6,f_ptj7,f_ml1j1,f_dphil1l2,f_maxmjoverpt,f_ptl1,f_detal1l2,f_ptj8,f_ptl3);
                f_nbtags = ss::bdt_jec_dn_nbtags();
                f_njets = ss::bdt_jec_dn_njets();
                f_met = ss::bdt_jec_dn_met();
                f_htb = ss::bdt_jec_dn_htb();
                f_nlb40 = ss::bdt_jec_dn_nlb40();
                f_ntb40 = ss::bdt_jec_dn_ntb40();
                mvavaluedn = get_prediction(f_nbtags,f_njets,f_met,f_ptl2,f_nlb40,f_ntb40,f_nleps,f_htb,f_q1,f_ptj1,f_ptj6,f_ptj7,f_ml1j1,f_dphil1l2,f_maxmjoverpt,f_ptl1,f_detal1l2,f_ptj8,f_ptl3);
                f_nbtags = ss::bdt_jer_up_nbtags();
                f_njets = ss::bdt_jer_up_njets();
                f_met = ss::bdt_jer_up_met();
                f_htb = ss::bdt_jer_up_htb();
                f_nlb40 = ss::bdt_jer_up_nlb40();
                f_ntb40 = ss::bdt_jer_up_ntb40();
                mvavalueJERup = get_prediction(f_nbtags,f_njets,f_met,f_ptl2,f_nlb40,f_ntb40,f_nleps,f_htb,f_q1,f_ptj1,f_ptj6,f_ptj7,f_ml1j1,f_dphil1l2,f_maxmjoverpt,f_ptl1,f_detal1l2,f_ptj8,f_ptl3);
                f_nbtags = ss::bdt_jer_dn_nbtags();
                f_njets = ss::bdt_jer_dn_njets();
                f_met = ss::bdt_jer_dn_met();
                f_htb = ss::bdt_jer_dn_htb();
                f_nlb40 = ss::bdt_jer_dn_nlb40();
                f_ntb40 = ss::bdt_jer_dn_ntb40();
                mvavalueJERdn = get_prediction(f_nbtags,f_njets,f_met,f_ptl2,f_nlb40,f_ntb40,f_nleps,f_htb,f_q1,f_ptj1,f_ptj6,f_ptj7,f_ml1j1,f_dphil1l2,f_maxmjoverpt,f_ptl1,f_detal1l2,f_ptj8,f_ptl3);
            }

            // float _event_bdt = numeric_limits<float>::quiet_NaN();
            // auto event_bdt = [&]() {
            //     if (!evaluateBDT) return 0.0f;
            //     if (isnan(_event_bdt)) _event_bdt = reader.EvaluateMVA("BDT");
            //     return _event_bdt;
            // };

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
                do_fill(h_nlb40, nlb40, (ss::is_real_data() ? 1 : 1/ss::weight_btagsf())); // Test
                do_fill(h_ntb40, ntb40);
                do_fill(h_nbtags, nbtags);
                do_fill(h_nbtagsheavyup, nbtags, ss::weight_btagsf_heavy_UP()/ss::weight_btagsf());
                do_fill(h_nbtagsheavydown, nbtags, ss::weight_btagsf_heavy_DN()/ss::weight_btagsf());
                do_fill(h_nbtagslightup, nbtags, ss::weight_btagsf_light_UP()/ss::weight_btagsf());
                do_fill(h_nbtagslightdown, nbtags, ss::weight_btagsf_light_DN()/ss::weight_btagsf());
                do_fill(h_bdisc1, bdisc1);
                do_fill(h_maxmjoverpt, maxmjoverpt);
                // if (ptj1>500) {
                //     do_fill(h_dphimetj1, calcDeltaPhi(jet1phi,metphi));
                //     do_fill(h_metoverptj1, met/ptj1);
                //     // calculate mht
                //     float ht_recalc = 0.;
                //     Vec4 mht_p4(0,0,0,0);
                //     for (int ijet = 0; ijet < ss::jets().size(); ijet++) {
                //         auto jet = ss::jets()[ijet]*ss::jets_undoJEC()[ijet]*ss::jets_JEC()[ijet];
                //         ht_recalc += jet.pt();
                //         mht_p4 -= jet;
                //     }
                //     mht_p4 -= ss::lep1_p4();
                //     mht_p4 -= ss::lep2_p4();
                //     TVector2 mhtVector = TVector2(mht_p4.pt()*cos(mht_p4.phi()), mht_p4.pt()*sin(mht_p4.phi()));
                //     TVector2 metVector = TVector2(ss::met()*cos(ss::metPhi()), ss::met()*sin(ss::metPhi()));
                //     float diffMetMht = (mhtVector - metVector).Mod();
                //     do_fill(h_metmmht,(diffMetMht/ss::met()));
                // }

                int nbnj = 5*min(nbtags,3)+(max(min(njets,6),2)-2);
                do_fill(h_nbnj, nbnj);

                int btagid = 0;
                // for (auto flav : ss::btags_flavor()) {
                //     // 20*b + 5*c + light
                //     if (abs(flav)==5) btagid += 20;
                //     else if (abs(flav)==4) btagid += 5;
                //     else btagid += 1;
                // }
                if (nbtags == 3) {
                    do_fill(h_btagid, btagid);
                }
                // if (ss::extragenb() == 2)

                do_fill(h_pt1, lep1ccpt);
                do_fill(h_pt2, lep2ccpt);
                if (nleps > 2) do_fill(h_pt3, lep3pt);
                do_fill(abs(lep1id) == 11 ? h_pte      : h_ptm,      lep1ccpt);
                do_fill(abs(lep2id) == 11 ? h_pte      : h_ptm,      lep2ccpt);

                do_fill(h_eta1,   ss::lep1_p4().eta());
                do_fill(h_eta2,   ss::lep2_p4().eta());
                do_fill(abs(lep1id) == 11 ? h_etae     : h_etam,     lep1eta);
                do_fill(abs(lep2id) == 11 ? h_etae     : h_etam,     lep2eta);
                do_fill(abs(lep1id) == 11 ? h_phie     : h_phim,     lep1phi);
                do_fill(abs(lep2id) == 11 ? h_phie     : h_phim,     lep2phi);
                do_fill(abs(lep1id) == 11 ? h_cphie     : h_cphim,     lep1phi);
                do_fill(abs(lep2id) == 11 ? h_cphie     : h_cphim,     lep2phi);

                if (abs(lep1id) == 11) do_fill(h_nmiss1, nmiss1);
                if (abs(lep2id) == 11) do_fill(h_nmiss2, nmiss2);

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

                bool lep1_fail_ptratio = false;
                bool lep2_fail_ptratio = false;

                if (is2017 || is2018) {
                    lep1_fail_ptratio = (lep1ptratio >= 0.) && (lep1ptratio < (abs(lep1id) == 11 ? 0.85 : 0.80));
                    lep2_fail_ptratio = (lep2ptratio >= 0.) && (lep2ptratio < (abs(lep2id) == 11 ? 0.85 : 0.80));
                } else {
                    lep1_fail_ptratio = (lep1ptratio >= 0.) && (lep1ptratio < (abs(lep1id) == 11 ? 0.80 : 0.76));
                    lep2_fail_ptratio = (lep2ptratio >= 0.) && (lep2ptratio < (abs(lep2id) == 11 ? 0.80 : 0.76));
                }

                if (lep1_fail_ptratio) do_fill(h_ptrelfail1, lep1ptrel);
                if (lep2_fail_ptratio) do_fill(h_ptrelfail2, lep2ptrel);
                if (lep1_fail_ptratio) do_fill(abs(lep1id) == 11 ? h_ptrelfaile   : h_ptrelfailm,   lep1ptrel);
                if (lep2_fail_ptratio) do_fill(abs(lep2id) == 11 ? h_ptrelfaile   : h_ptrelfailm,   lep2ptrel);
                if (lep1_fail_ptratio and looseleg == 1) do_fill(h_ptrelfaillnt, lep1ptrel);
                if (lep2_fail_ptratio and looseleg == 2) do_fill(h_ptrelfaillnt, lep2ptrel);

                do_fill(h_ptratio1, lep1ptratio);
                do_fill(h_ptratio2, lep2ptratio);
                if (looseleg > 0) do_fill(h_ptratiolnt, looseleg == 1 ? lep1ptratio   : lep2ptratio);
                do_fill(abs(lep1id) == 11 ? h_ptratioe : h_ptratiom, lep1ptratio);
                do_fill(abs(lep2id) == 11 ? h_ptratioe : h_ptratiom, lep2ptratio);

                do_fill(h_miniiso1, lep1miniiso);
                do_fill(h_miniiso2, lep2miniiso);
                if (looseleg > 0) do_fill(h_miniisolnt, looseleg == 1 ? lep1miniiso   : lep2miniiso);
                do_fill(abs(lep1id) == 11 ? h_miniisoe : h_miniisom, lep1miniiso);
                do_fill(abs(lep2id) == 11 ? h_miniisoe : h_miniisom, lep2miniiso);

                do_fill(h_dphil1l2, ss::bdt_dphil1l2());
                do_fill(h_detal1l2, ss::bdt_detal1l2());
                do_fill(h_absdetal1l2, ss::bdt_detal1l2());


                do_fill(h_matchtype,   matchtype);

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

                if (njets >= 1) do_fill(h_ptj1, ptj1);
                if (njets >= 2) do_fill(h_ptj2, ptj2);
                if (njets >= 3) do_fill(h_ptj3, ptj3);
                if (njets >= 4) do_fill(h_ptj4, ptj4);
                if (not doSS) {
                    if (njets >= 5) do_fill(h_ptj5, ptj5);
                    if (njets >= 6) do_fill(h_ptj6, ptj6);
                    if (njets >= 7) do_fill(h_ptj7, ptj7);
                    if (njets >= 8) do_fill(h_ptj8, ptj8);

                    if (njets >= 1) do_fill(h_ml1j1, ml1j1);



                    // // do_fill(h_event_bdt, event_bdt());
                    // do_fill(h_event_bdt, ss::bdt_disc());
                    // do_fill(h_event_bdt_scaleup, ss::bdt_disc(), doScaleUnc ? 1.0 : ss::weight_scale_UP());
                    // do_fill(h_event_bdt_scaledn, ss::bdt_disc(), doScaleUnc ? 1.0 : ss::weight_scale_DN());
                    // do_fill(h_event_bdt_btagheavyup, ss::bdt_disc(), ss::is_real_data() ? 1.0 : ss::weight_btagsf_heavy_UP()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_btagheavydn, ss::bdt_disc(), ss::is_real_data() ? 1.0 : ss::weight_btagsf_heavy_DN()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_btaglightup, ss::bdt_disc(), ss::is_real_data() ? 1.0 : ss::weight_btagsf_light_UP()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_btaglightdn, ss::bdt_disc(), ss::is_real_data() ? 1.0 : ss::weight_btagsf_light_DN()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_jecup, ss::bdt_disc_jec_up());
                    // do_fill(h_event_bdt_jecdn, ss::bdt_disc_jec_dn());
                    // do_fill(h_event_bdt_jerup, ss::bdt_disc_jer_up());
                    // do_fill(h_event_bdt_jerdn, ss::bdt_disc_jer_dn());
                    //
                    // do_fill(h_event_bdt, event_bdt());
                    // do_fill(h_event_bdt, mvavalue);
                    // do_fill(h_event_bdt_scaleup, mvavalue, doScaleUnc ? 1.0 : ss::weight_scale_UP());
                    // do_fill(h_event_bdt_scaledn, mvavalue, doScaleUnc ? 1.0 : ss::weight_scale_DN());
                    // do_fill(h_event_bdt_btagheavyup, mvavalue, ss::is_real_data() ? 1.0 : ss::weight_btagsf_heavy_UP()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_btagheavydn, mvavalue, ss::is_real_data() ? 1.0 : ss::weight_btagsf_heavy_DN()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_btaglightup, mvavalue, ss::is_real_data() ? 1.0 : ss::weight_btagsf_light_UP()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_btaglightdn, mvavalue, ss::is_real_data() ? 1.0 : ss::weight_btagsf_light_DN()/ss::weight_btagsf());
                    // do_fill(h_event_bdt_jecup, mvavalueup);
                    // do_fill(h_event_bdt_jecdn, mvavaluedn);
                    // do_fill(h_event_bdt_jerup, mvavalueJERup);
                    // do_fill(h_event_bdt_jerdn, mvavalueJERdn);

                    // float adhoc = 1.;
                    // if (nrealb >= 3) adhoc = 1.35;
                    // do_fill(h_event_bdt_bbup, mvavalue, ss::is_real_data() ? 1 : adhoc);
                }
            };

            bool BR_lite = ht > 300       and njets >= 2 and
                           nbtags >= 2    and met >= 50;
            bool BR = BR_lite and hyp_class == 3;
            bool CRW = BR and nleps == 2 and nbtags == 2 and njets <= 5;
            bool CRZ = BR_lite and hyp_class == 6 and (zcand13 or zcand23) and lep3ccpt > 20 and
                       (class6Fake or (lep1good and lep2good and lep3good));

            if (CRZ) fill_region("crz", weight);

            if (CRW) fill_region("crw", weight);

            if (BR) fill_region("br", weight);

            if (hyp_class == 4) fill_region("osloose", weight);
            if (hyp_class == 4 and lep1ccpt > 25 and lep2ccpt > 25) fill_region("hhosloose", weight);

            if (hyp_class == 3 and njets >= 2) {
                fill_region("nj2", weight);
            }

            if (hyp_class == 3 and nbtags == 1 and
                   njets >= 2 and met > 50) {
                fill_region("htnb1", weight);

                // if ( (abs(lep1id) == 13 and (lep1pt < 150 and lep1pt > 125)) or
                //      (abs(lep2id) == 13 and (lep2pt < 150 and lep2pt > 125))
                //    ) {
                //     // if (ss::is_real_data() and !doFakes and !doFlips) {
                //     //     cout << ss::run() << ":" << ss::lumi() << ":" << ss::event() << " " << lep1id << " " << lep2id << " " << lep1pt << " " << lep2pt << " " << ss::filename() << endl;
                //     // }
                //     fill_region("htnb1mu", weight);
                // }

            }

            if (hyp_class == 3 and njets >= 2 and met > 50
                    and lep1ccpt > 25 and lep2ccpt > 25
                    ) {
                fill_region("sshh", weight);
            }

            if (njets >= 2 and
                    met > 50. and
                    ht > 300) {
                if (hyp_class == 4 and !zcand12) {
                    fill_region("os_noisr", weight / isr_reweight(useIsrWeight, year, nisrmatch));
                    fill_region("os", weight);
                    // if (lep1ccpt > 25 and lep2ccpt > 20) {
                    //     fill_region("ftos", weight);
                    // }
                    fill_region("osnbrw", weight * nb_reweight(nbtags));
                    if (lep1ccpt > 25 and lep2ccpt > 25) fill_region("hhos", weight);
                    // if (evaluateBDT) {
                    //     if (event_bdt() > 0.21) {
                    //         fill_region("os_highbdt", weight);
                    //     } else {
                    //         fill_region("os_lowbdt", weight);
                    //     }
                    // }
                }
                if (hyp_class == 2) fill_region("tl", weight);
                if (hyp_class == 2 and lep1ccpt > 25 and lep2ccpt > 25) fill_region("hhtl", weight);
            }

            if (njets >= 2 and
                    met > 50. and
                    hyp_class == 4 and !zcand12) {
                fill_region("os_noht", weight);
            }

            if (njets >= 2 and
                    met > 50. and
                    ht > 300) {
                if (hyp_class == 4 and !zcand12 and nbtags>=3) {
                    fill_region("oshtnb3", weight);
                }
            }


            if (hyp_class == 3 and m12 < 76 and njets >= 2 and nbtags == 0 and nleps >= 3) {
                if (m3l > 50. and m3l < 120.) {
                    fill_region("xgcr", weight);
                }
            }

            // for SS
            if (lep1ccpt > 25. and lep2ccpt > 25.) {
                if (not (m12 > 76 and m12 < 106 and abs(lep1id)==11 and abs(lep2id)==11)) {
                    if (njets >= 2 and met < 50. and ht < 300. and hyp_class == 3) fill_region("lowmetlowht", weight);
                    if (njets >= 2 and met < 50.               and hyp_class == 3) fill_region("lowmetallht", weight);
                }
            }
            if (njets >= 2 and               hyp_class == 4) fill_region("osnomet", weight);
            if (njets >= 2 and               hyp_class == 4 and ss::mtmin() > 150.) fill_region("osnomethighmt", weight);
            if (njets >= 2 and               hyp_class == 2) fill_region("tlnomet", weight);
            if (njets >= 2 and met > 50. and hyp_class == 4) fill_region("osmet", weight);
            if (njets >= 2 and met > 500. and hyp_class == 4) fill_region("osmethigh", weight);
            // if (njets >= 2 and met > 200. and hyp_class == 4 and met/ptj1>0.3 and ptj1>500) {
            //     fill_region("osmetjet", weight);
            // }
            if (njets >= 2 and met > 50. and hyp_class == 2) fill_region("tlmet", weight);
            if (njets >= 2 and met < 50. and hyp_class == 6 and // NOTE LT
                    lep1ccpt > 25. and
                    ((lep1good and lep2good and lep3good) or class6Fake)
                    ) {
                // if (lep2ccpt > 20. and lep3ccpt > 20.) {
                // }

                // for 3 leptons, trigger is pretty much fully efficient, so divide out the scale factor we multiplied in earlier
                fill_region("mlonz", weight/triggerScaleFactor(year, lep1id, lep2id, lep1pt, lep2pt, lep1eta, lep2eta, ht, analysis, 0));
                // fill_region("mlonz", weight);
            }
            // if (njets >= 2 and met > 50. and hyp_class == 6 and // NOTE GT
            //         lep1ccpt > 25. and
            //         lep1good and lep2good and lep3good
            //    ) {
            //     fill_region("mlonzhigh", weight);
            // }
            // if (njets >= 0 and met > 50. and hyp_class == 6 and // NOTE GT
            //         lep1ccpt > 25. and
            //         lep1good and lep2good and lep3good
            //    ) {
            //     fill_region("mlonz0j", weight);
            // }

            // // BDT Validation Regions
            // if (hyp_class == 3 and
            //         njets >= 2) {
            //     if (nbtags < 2  and ht > 300 and met > 50) fill_region("bdt_nb", weight);  // invert Nb
            //     if (nbtags >= 2 and ht < 300 and met > 50) {                               // invert Ht
            //         fill_region("bdt_ht", weight);
            //         fill_region("bdt_met_ht", weight);
            //     }
            //     if (nbtags >= 2 and ht > 300 and met < 50) {                               // invert MET
            //         fill_region("bdt_met", weight);
            //         fill_region("bdt_met_ht", weight);
            //     }
            // }

            // bool BDT_train = hyp_class == 3 and
            //                  ht > 250       and njets >= 2    and
            //                  nbtags >= 1;
            // if ((BDT_train and not BR) or CRW) fill_region("bdt_train", weight);

            if (hyp_class == 4 and !zcand12 and (njets-nisrjets) == 2 and
                    lep1ccpt > 25 and lep2ccpt > 20 and njets >= 2) {
                fill_region("tt_isr", weight / isr_reweight(useIsrWeight, year, nisrmatch));
                fill_region("tt_isr_reweight_check", weight);
            }
            // if (hyp_class == 4 and !zcand12 and (njets-nisrjets) == 2 and
            //         lep1ccpt > 25 and lep2ccpt > 20 and njets >= 2) {
            //     fill_region("tt_isr_test", weight / isr_reweight(useIsrWeight, year, nisrmatch));
            // }

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

