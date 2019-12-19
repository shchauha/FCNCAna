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
#include "TLorentzVector.h"
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
  map<string, TH1D> ee;
  map<string, TH1D> em;
  map<string, TH1D> mm;

  HistCol(vector<string> regions, const string& name, int nbins, const float* bins, vector<HistCol*>* registry=nullptr) {
    for (string region : regions) {
      string base_name = region + "_" + name;
      string base_title = region + " " + name;
      in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, bins));
      ee.emplace(region, TH1D((base_name + "_ee").c_str(), (base_title + " ee").c_str(), nbins, bins));
      em.emplace(region, TH1D((base_name + "_em").c_str(), (base_title + " em").c_str(), nbins, bins));
      mm.emplace(region, TH1D((base_name + "_mm").c_str(), (base_title + " mm").c_str(), nbins, bins));

    }
    if (registry != nullptr)
      registry->push_back(this);
  }

  HistCol(vector<string> regions, const string& name, int nbins, float low, float high, vector<HistCol*>* registry=nullptr) {
    for (string region : regions) {
      string base_name = region + "_" + name;
      string base_title = region + " " + name;
      in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, low, high));
      ee.emplace(region, TH1D((base_name + "_ee").c_str(), (base_title + " ee").c_str(), nbins, low, high));
      em.emplace(region, TH1D((base_name + "_em").c_str(), (base_title + " em").c_str(), nbins, low, high));
      mm.emplace(region, TH1D((base_name + "_mm").c_str(), (base_title + " mm").c_str(), nbins, low, high));
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
    for (auto p : in) p.second.Write();
    for (auto p : ee) p.second.Write();
    for (auto p : em) p.second.Write();
    for (auto p : mm) p.second.Write();

  }
};

int signal_region(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps) {
  // remember that sgn(pdgid) != sgn(charge), it's flipped. so mad.
  int mm = (id1 > 0);
  if ( met>300 && ht<300) return -1;
  //High-high
  if (met >= 300 && ht >= 300) {
    if (njets <= 4) {
      if (met < 500) return 46+mm;
      else return 48+mm;
    } else {
      if (met < 500) return 50+mm;
      else return 52+mm;
    }
  }
  if (ht >= 1125) {
    if (ht < 1300) {
      if (njets <= 4) return 54;
      else if (njets <=6) return 57;
      else return 60;
    }
    else if (ht < 1600) {
      if (njets <= 4) return 55;
      else if (njets <=6) return 58;
      else return 61;
    }
    else {
      if (njets <= 4) return 56;
      else if (njets <=6) return 59;
      else return 62;
    }
  }
  if (ht < 300) {
    if (nbtags == 0 && mt_min < 120 && met < 200 && njets <= 4) return 1; 
    if (nbtags == 0)return 3;     
    if (nbtags == 1 && mt_min < 120 && met < 200 && njets <= 4) return 11;
    if (nbtags == 1) return 13+mm; 
    if (nbtags == 2 && mt_min < 120 && met < 200 && njets <= 4) return 23; 
    if (nbtags == 2) return 25+mm; 
    if (nbtags >= 3) {
      if (mt_min < 120) return 35+mm; 
      else return 41;
    }
  }
  if (ht >= 300 && ht < 1125) {
    if (nbtags == 0){
      if (mt_min < 120 && met < 200 && njets <= 4) return 2; 
      if (mt_min < 120 && met < 200 && njets > 4) return 4; 
      if (mt_min < 120 && met >= 200 && njets <= 4) return 5+mm; 
      if (mt_min < 120 && met >= 200 && njets > 4) return 7; 
      if (mt_min >= 120 && met < 200 && njets <= 4) return 8+mm;
      return 10;
    } 
    if (nbtags == 1){
      if (mt_min < 120 && met < 200 && njets <= 4) return 12; 
      if (mt_min < 120 && met < 200 && njets > 4) return 15+mm; 
      if (mt_min < 120 && met >= 200 && njets <= 4) return 17+mm; 
      if (mt_min < 120 && met >= 200 && njets > 4) return 19; 
      if (mt_min >= 120 && met < 200 && njets <= 4) return 20+mm;
      return 22;
    } 
    if (nbtags == 2){
      if (mt_min < 120 && met < 200 && njets <= 4) return 24; 
      if (mt_min < 120 && met < 200 && njets > 4) return 27+mm; 
      if (mt_min < 120 && met >= 200 && njets <= 4) return 29+mm; 
      if (mt_min < 120 && met >= 200 && njets > 4) return 31; 
      if (mt_min >= 120 && met < 200 && njets <= 4) return 32+mm;
      return 34;
    } 
    if (nbtags >= 3){
      if (mt_min < 120) {
	if (njets <= 4) return 37+mm;
	else return 39+mm;
      } else {
	if (njets <= 4) return 42+mm;
	else return 44+mm;
      }
    }
  }
 
  return -1;
 
}

float calcDeltaPhi(float phi1, float phi2){
  float dPhi = phi1 - phi2;
  while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
  while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
  return fabs(dPhi);
}

float calcDeltaR(float eta1, float phi1, float eta2, float phi2){
  return TMath::Sqrt(pow(calcDeltaPhi(phi1, phi2),2)+pow((eta1-eta2),2));
  
}

float calcMT(float pt1, float phi1, float pt2, float phi2){
  return sqrt( 2 * pt1 * pt2 * ( 1 - cos( phi1 - phi2 ) ) );
}


float minDR(const Vec4 & lep , const vector<Vec4> & jet)
{
  int size = (int)jet.size();
  //cout<<"jet size "<<size<<endl;
  float mindr = 999;  
  if(size)for(int i=0;i<size;i++){
      float dr = calcDeltaR(lep.eta(), lep.phi(), jet[i].eta(), jet[i].phi());
      if(dr<mindr) mindr = dr;
    }
  return mindr;
}

float PtMaxEta(const vector<Vec4> & jet)
{
  int size = (int)jet.size();
  //cout<<"jet size "<<size<<endl;
  float maxeta=0, eta=0;
  int index=0;
  for(int i=0;i<size;i++){
        eta = jet[i].eta();
        if(eta>maxeta){
                maxeta = eta;
                index = i;
        }
  }
  return jet[index].pt();
}

struct lepton
{
  Vec4 v;
  int id;
  bool isgood;
  float miniiso;
  float dxy;
  float dz;   
  
}; 

float MOSSF(vector<lepton> lep)
{
  int size = (int)lep.size();
  float mossf = 999999;
  float diff = 999999;
  if(size>=2)for(int i=0;i<size;i++){
      for(int j=i+1;j<size;j++){	
	if(lep[i].id*lep[j].id==-169||lep[i].id*lep[j].id==-121){
	  float mass = (lep[i].v+lep[j].v).M();	  
	  if(TMath::Abs(mass-91.2)<diff)
	    {
	      mossf = mass;
	      diff = TMath::Abs(mass-91.2);
	  }
	}	  
      }
    }
  
  return mossf;
}


float isr_reweight(bool useIsrWeight, int year, int nisrmatch, int sample) {
  if (!useIsrWeight) return 1;
  if (ss::is_real_data()) return 1;
  return isrWeight(year, nisrmatch, sample); // 10 is ttbar
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
  bool BDTTraining = options.Contains("BDTTraining");
  bool BDTApplication = options.Contains("BDTApplication");
  bool ReadBDT = options.Contains("ReadBDT");
  bool isData = options.Contains("doData");
  
  //ana_t analysis = FTANA;
  //if (doSS) {
  ana_t analysis = SSANA;
  //}

  TString proc(ch->GetTitle());

  //cout<<"proc "<<proc<<endl;
  // bool useIsrWeight = proc.Contains("tt_");
  bool useIsrWeight = proc.Contains("tt_") 
    or proc.Contains("ttw") 
    or proc.Contains("ttz") 
    or proc.Contains("tth");


  bool doScaleUnc = proc.Contains("tt_");
  bool useTTBB = proc.Contains("tt_") 
    or proc.Contains("ttw") 
    or proc.Contains("ttz") 
    or proc.Contains("tth");
  
  useTTBB = false;
  //cout<<"using isr "<<useIsrWeight<<" TTBB "<<useTTBB<<endl;  
  if (noISRWeights) useIsrWeight = false;
  // We may have derived the fake rate map throwing away leptons with pT<18 (e.g., 2017), so
  // we need to apply this cut here to be consistent
  //float min_pt_fake = minPtFake18 ? 18. : -1;
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

  vector<string> regions =
    {			    
      //"sshh", 
      "ssbr",
      "ss0b2j",
      "ss1b2j",
      "ss2b2j",
      "mlbr",
      "ml1b1j",
      "ml2b2j",
      // "mlbronz",
      // "ml1b1jonz",
      // "ml2b2jonz",
      // "mlbrinc",
      // "ml1b1jinc",
      // "ml2b2jinc",
      "osbr",
      "tl",
      
      //"br",
      //"susytl",
      //"hhtl",     
      //"lowmetonzor0b",
      //"mllowmetonz2b",
      
      // "ssbr2",
      // "ss1b2j2",
      // "ss2b2j2",
      // "ss1b2jbtagM",
      // "ss2b2jbtagM",
      // "ss1b2jbtag25",
      // "ss2b2jbtag25",
      // "ss1b2jbtag25M",
      // "ss2b2jbtag25M",
      // "ss1b2jjet40",
      // "ss2b2jjet40",
      // "ss1b2jjet40btagM",
      // "ss2b2jjet40btagM",
      // "ss1b2jjet40btag25",
      // "ss2b2jjet40btag25",
      // "ss1b2jjet40btag25M",
      // "ss2b2jjet40btag25M",    
      // //"ssonz",
      // //"ss1b2j_s",
      // //"ss1b2j_s_met50",
    
      
    };

  vector<HistCol*> registry;
  vector<HistCol2D*> registry2D;
  // HistCol h_met         (regions, "met"        , 30, 0   , 300 , &registry);
  HistCol h_met           (regions, "met"             , 20, 0   , 500 ,  &registry);
  HistCol h_metphi        (regions, "metphi"          , 35, -3.5, 3.5 ,  &registry);
  HistCol h_rawmet        (regions, "rawmet"          , 20, 0   , 500 ,  &registry);
  //HistCol h_calomet      (regions, "calomet"         , 60, 0   , 600 ,  &registry);
  HistCol h_ht            (regions, "ht"              , 20, 0   , 1000,  &registry);
  HistCol h_mll           (regions, "mll"             , 30, 0   , 300 ,  &registry);
  HistCol h_mllzoom       (regions, "mllzoom"         , 12, 60  , 120 ,  &registry);
  HistCol h_m3l           (regions, "m3l"             , 20, 0   , 200 ,  &registry);
  HistCol h_mtmin         (regions, "mtmin"           , 15, 0   , 300 ,  &registry);
  HistCol h_dphil1met     (regions, "dphil1met"       , 35, -3.5, 3.5 ,  &registry);
  HistCol h_dphil2met     (regions, "dphil2met"       , 35, -3.5, 3.5 ,  &registry);
  HistCol h_dphimetj1     (regions, "dphimetj1"       , 35, -3.5, 3.5 ,  &registry);
  HistCol h_metoverptj1   (regions, "metoverptj1"     , 50, 0.,2.,       &registry);
  HistCol h_nleps         (regions, "nleps"           , 5, -0.5 , 4.5 ,  &registry);
  HistCol h_nele          (regions, "nele"            , 5, -0.5 , 4.5 ,  &registry);
  HistCol h_njets         (regions, "njets"           , 6 , 0   , 6   ,  &registry);
  HistCol h_njets40       (regions, "njets40"         , 6 , 0   , 6   ,  &registry);
  HistCol h_nisrjets      (regions, "nisrjets"        , 5 , 0   , 5   ,  &registry);
  HistCol h_nisrmatch     (regions, "nisrmatch"       , 5 , 0   , 5   ,  &registry);
  HistCol h_nlb40         (regions, "nlb40"           , 5 , 0   , 5   ,  &registry);
  HistCol h_ntb40         (regions, "ntb40"           , 8 , 0   , 8   ,  &registry);
  HistCol h_nbtags        (regions, "nbtags"          , 5 , 0   , 5   ,  &registry);
  HistCol h_nbtagszoom    (regions, "nbtagszoom"      , 3 , 0   , 3   ,  &registry);
  HistCol h_nbtagsM       (regions, "nbtagsM"         , 5 , 0   , 5   ,  &registry);
  HistCol h_nbtags25      (regions, "nbtags25"        , 5 , 0   , 5   ,  &registry);
  HistCol h_nbtags25M     (regions, "nbtags25M"       , 5 , 0   , 5   ,  &registry);
  HistCol h_maxmjoverpt   (regions, "maxmjoverpt"     , 50, 0   , 0.35,  &registry);
  HistCol h_btagid        (regions, "btagid"          , 100,0   , 100 ,  &registry);
  HistCol h_pt1           (regions, "pt1"             , 30, 0   , 300 ,  &registry);
  HistCol h_pt2           (regions, "pt2"             , 30, 0   , 300 ,  &registry);
  HistCol h_pt3           (regions, "pt3"             , 30, 0   , 300 ,  &registry);
  HistCol h_pte           (regions, "pte"             , 30, 0   , 300 ,  &registry);
  HistCol h_ptm           (regions, "ptm"             , 30, 0   , 300 ,  &registry);
  HistCol h_eta1          (regions, "eta1"            , 35, -3.5, 3.5 ,  &registry);
  HistCol h_eta2          (regions, "eta2"            , 35, -3.5, 3.5 ,  &registry);
  HistCol h_etae          (regions, "etae"            , 35, -3.5, 3.5 ,  &registry);
  HistCol h_etam          (regions, "etam"            , 35, -3.5, 3.5 ,  &registry);

  HistCol h_ptrel1        (regions, "ptrel1"          , 15, 0   , 50  ,  &registry);
  HistCol h_ptrel2        (regions, "ptrel2"          , 15, 0   , 50  ,  &registry);
  HistCol h_ptrele        (regions, "ptrele"          , 15, 0   , 50  ,  &registry);
  HistCol h_ptrelm        (regions, "ptrelm"          , 15, 0   , 50  ,  &registry);

  HistCol h_ptratio1      (regions, "ptratio1"        , 30, 0   , 1.5 ,  &registry);
  HistCol h_ptratio2      (regions, "ptratio2"        , 30, 0   , 1.5 ,  &registry);
  HistCol h_ptratioe      (regions, "ptratioe"        , 30, 0   , 1.5 ,  &registry);
  HistCol h_ptratiom      (regions, "ptratiom"        , 30, 0   , 1.5 ,  &registry);

  HistCol h_miniiso1      (regions, "miniiso1"        , 15,  0  , 0.2 ,  &registry);
  HistCol h_miniiso2      (regions, "miniiso2"        , 15,  0  , 0.2 ,  &registry);
  HistCol h_miniisoe      (regions, "miniisoe"        , 15,  0  , 0.2 ,  &registry);
  HistCol h_miniisom      (regions, "miniisom"        , 15,  0  , 0.2 ,  &registry);

  HistCol h_dphil1l2      (regions, "dphil1l2"        , 15,  0  , 4   ,  &registry);
  HistCol h_detal1l2      (regions, "detal1l2"        , 30,  -4 , 4   ,  &registry);
  HistCol h_absdetal1l2   (regions, "absdetal1l2"     , 15,  0  , 4   ,  &registry);

  HistCol h_type          (regions, "type"            , 3 , -0.5, 2.5 ,  &registry);
  HistCol h_type3l        (regions, "type3l"          , 4 , -0.5, 3.5 ,  &registry);
  HistCol h_nvtx          (regions, "nvtx"            , 31, -0.5, 61.5,  &registry);

  HistCol h_ptj1          (regions, "ptj1"            , 50, 0   , 500 ,  &registry);
  HistCol h_ptj2          (regions, "ptj2"            , 50, 0   , 500 ,  &registry);
  HistCol h_ptj3          (regions, "ptj3"            , 50, 0   , 500 ,  &registry);
  HistCol h_ptj4          (regions, "ptj4"            , 50, 0   , 500 ,  &registry);
  HistCol h_ptj5          (regions, "ptj5"            , 50, 0   , 500 ,  &registry);
  HistCol h_ptj6          (regions, "ptj6"            , 30, 0   , 300 ,  &registry);
  HistCol h_ptj7          (regions, "ptj7"            , 30, 0   , 300 ,  &registry);
  HistCol h_ptj8          (regions, "ptj8"            , 30, 0   , 300 ,  &registry);
  
  HistCol h_ptbt1         (regions, "ptbt1"           , 50, 0   , 500 ,  &registry);
  HistCol h_ptbt2         (regions, "ptbt2"           , 50, 0   , 500 ,  &registry);
  HistCol h_ptbt3         (regions, "ptbt3"           , 50, 0   , 500 ,  &registry);
  HistCol h_ptbt4         (regions, "ptbt4"           , 50, 0   , 500 ,  &registry);
  HistCol h_nforwardjets20(regions, "nforwardjets20"  , 10, 0 , 10 ,     &registry);
  HistCol h_sr            (regions, "sr"              , 60, -0.5, 59.5,  &registry);
  HistCol h_bdt           (regions, "bdt"             , 10 , -1., 1.,    &registry);
  HistCol h_bdt_hut       (regions, "bdt_hut"         , 10 , -1., 1.,    &registry);
  HistCol h_bdt_hut_ttbar (regions, "bdt_hut_ttbar"   , 10 , -1., 1.,    &registry);
  HistCol h_bdt_hut_ttv   (regions, "bdt_hut_ttv"     , 10 , -1., 1.,    &registry);
  HistCol h_bdt_hct       (regions, "bdt_hct"         , 10 , -1., 1.,    &registry);
  HistCol h_bdt_hct_ttbar (regions, "bdt_hct_ttbar"   , 10 , -1., 1.,    &registry);
  HistCol h_bdt_hct_ttv   (regions, "bdt_hct_ttv"     , 10 , -1., 1.,    &registry);
  
  HistCol h_drl1l2        (regions, "drl1l2"          , 20,  0  , 5   ,  &registry);
  HistCol h_mindrl1j      (regions, "mindrl1j"        , 20,  0  , 5   ,  &registry);
  HistCol h_mindrl2j      (regions, "mindrl2j"        , 20,  0  , 5   ,  &registry);
  HistCol h_mindrl1bt     (regions, "mindrl1bt"       , 20,  0  , 5   ,  &registry);
  HistCol h_mindrl2bt     (regions, "mindrl2bt"       , 20,  0  , 5   ,  &registry);
  HistCol h_mt1           (regions, "mt1"             , 30, 0   , 300 ,  &registry);
  HistCol h_mt2           (regions, "mt2"             , 30, 0   , 300 ,  &registry);
  HistCol h_l1dxy         (regions, "l1dxy"           , 20, -0.05, 0.05 ,&registry);
  HistCol h_l2dxy         (regions, "l2dxy"           , 20, -0.05, 0.05 ,&registry);
  HistCol h_l1dz          (regions, "l1dz"            , 20, -0.1, 0.1 ,  &registry);
  HistCol h_l2dz          (regions, "l2dz"            , 20, -0.1, 0.1 ,  &registry);
  HistCol h_mossf         (regions, "mossf"           , 30, 0, 300 ,     &registry);
  HistCol h_fwd_jetpt     (regions, "fwd_jetpt"       , 50, 0, 500 ,     &registry);


  // 2D hists
  HistCol2D h_bdt2D_hut         (regions, "bdt2D_hut"             , 5 , -1., 1., 5 , -1., 1.,    &registry2D);
  HistCol2D h_bdt2D_hct         (regions, "bdt2D_hct"             , 5 , -1., 1., 5 , -1., 1.,    &registry2D);

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

  int nleps, njets, njets40, nbtags, nbtagsM, nbtags25, nbtags25M ;
  float ht, htb, met, metphi, rawmet, calomet;

  float maxmjoverpt, ml1j1;
  int matchtype;
  int nlb40, ntb40, nisrjets, nisrmatch;
  float ptj1, ptj2, ptj3, ptj4, ptj5, ptj6, ptj7, ptj8;
  float ptbt1, ptbt2, ptbt3, ptbt4;
  float jet1phi;

  int SR;
  float weight;    
  int nEventsTotal = 0;
  int nEventsChain = ch->GetEntries();

  //add list of trees
  float tree_lep1pt = -1;
  float tree_lep2pt = -1;
  float tree_lep1eta = -1;
  float tree_lep2eta = -1;
  float tree_njets = -1;
  float tree_nbtags = -1;
  float tree_met = -1;
  float tree_ht = -1;
  float tree_mll = -1;
  float tree_nele = -1;
  float tree_jet1pt = -1;
  float tree_jet2pt = -1;
  float tree_btag1pt = -1;
  float tree_drl1l2 = -1;
  float tree_mindrl1j = -1;
  float tree_mindrl2j = -1;
  float tree_mindrl1bt = -1;
  float tree_mindrl2bt = -1;
  float tree_dphil1met = -1;
  float tree_dphil2met = -1;
  float tree_mt1 = -1;
  float tree_mt2 = -1;
  float tree_dphil1l2 = -1;
  float tree_l1miniiso = -1;
  float tree_l2miniiso = -1;
  float tree_l1dxy = -1;
  float tree_l1dz = -1;
  float tree_l2dxy = -1;
  float tree_l2dz = -1;
  float tree_l1ptratio = -1;
  float tree_l1ptrel = -1;
  float tree_l2ptratio = -1;
  float tree_l2ptrel = -1;
  float tree_weight = -1;
  float tree_jet3pt = -1;
  float tree_fwd_jetpt = -1;


  TFile *  f1 = new TFile(Form("%s/histos_%s.root", outputdir.Data(), ch->GetTitle()), "RECREATE");
  f1->cd();
  TTree* out_tree = new TTree("t","fortraining");

  out_tree->Branch("lep1pt", &tree_lep1pt );
  out_tree->Branch("lep2pt", &tree_lep2pt );
  out_tree->Branch("lep1eta", &tree_lep1eta );
  out_tree->Branch("lep2eta", &tree_lep2eta );
  out_tree->Branch("njets", &tree_njets );
  out_tree->Branch("nbtags", &tree_nbtags );
  out_tree->Branch("met", &tree_met );
  out_tree->Branch("ht", &tree_ht );
  out_tree->Branch("mll", &tree_mll );
  out_tree->Branch("nele", &tree_nele );
  out_tree->Branch("jet1pt", &tree_jet1pt );
  out_tree->Branch("jet2pt", &tree_jet2pt );
  out_tree->Branch("btag1pt", &tree_btag1pt );
  out_tree->Branch("drl1l2", &tree_drl1l2 );
  out_tree->Branch("mindrl1j", &tree_mindrl1j );
  out_tree->Branch("mindrl2j", &tree_mindrl2j );
  out_tree->Branch("mindrl1bt", &tree_mindrl1bt );
  out_tree->Branch("mindrl2bt", &tree_mindrl2bt );
  out_tree->Branch("dphil1met", &tree_dphil1met );
  out_tree->Branch("dphil2met", &tree_dphil2met );
  out_tree->Branch("mt1", &tree_mt1 );
  out_tree->Branch("mt2", &tree_mt2 );
  out_tree->Branch("dphil1l2", &tree_dphil1l2 );
  out_tree->Branch("l1miniiso", &tree_l1miniiso );
  out_tree->Branch("l2miniiso", &tree_l2miniiso );
  out_tree->Branch("l1dxy", &tree_l1dxy );
  out_tree->Branch("l1dz", &tree_l1dz );
  out_tree->Branch("l2dxy", &tree_l2dxy );
  out_tree->Branch("l2dz", &tree_l2dz );
  out_tree->Branch("l1ptratio", &tree_l1ptratio );
  out_tree->Branch("l1ptrel", &tree_l1ptrel );
  out_tree->Branch("l2ptratio", &tree_l2ptratio );
  out_tree->Branch("l2ptrel", &tree_l2ptrel );
  out_tree->Branch("jet3pt", &tree_jet3pt );
  out_tree->Branch("fwd_jetpt", &tree_fwd_jetpt );
  out_tree->Branch("weight", &tree_weight);

  
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  reader->AddVariable("lep1pt", &tree_lep1pt );
  reader->AddVariable("lep2pt", &tree_lep2pt );
  reader->AddVariable("lep1eta", &tree_lep1eta );
  reader->AddVariable("lep2eta", &tree_lep2eta );
  reader->AddVariable("njets", &tree_njets );
  reader->AddVariable("nbtags", &tree_nbtags );
  reader->AddVariable("met", &tree_met );
  reader->AddVariable("ht", &tree_ht );
  reader->AddVariable("mll", &tree_mll );
  reader->AddVariable("nele", &tree_nele );
  reader->AddVariable("jet1pt", &tree_jet1pt );
  reader->AddVariable("jet2pt", &tree_jet2pt );
  reader->AddVariable("btag1pt", &tree_btag1pt );
  reader->AddVariable("drl1l2", &tree_drl1l2 );
  reader->AddVariable("mindrl1j", &tree_mindrl1j );
  reader->AddVariable("mindrl2j", &tree_mindrl2j );
  reader->AddVariable("mindrl1bt", &tree_mindrl1bt );
  reader->AddVariable("mindrl2bt", &tree_mindrl2bt );
  reader->AddVariable("dphil1met", &tree_dphil1met );
  reader->AddVariable("dphil2met", &tree_dphil2met );
  reader->AddVariable("mt1", &tree_mt1 );
  reader->AddVariable("mt2", &tree_mt2 );
  reader->AddVariable("dphil1l2", &tree_dphil1l2 );
  reader->AddVariable("l1miniiso", &tree_l1miniiso );
  reader->AddVariable("l2miniiso", &tree_l2miniiso );
  reader->AddVariable("l1dxy", &tree_l1dxy );
  reader->AddVariable("l1dz", &tree_l1dz );
  reader->AddVariable("l2dxy", &tree_l2dxy );
  reader->AddVariable("l2dz", &tree_l2dz );
  reader->AddVariable("l1ptratio", &tree_l1ptratio );
  reader->AddVariable("l1ptrel", &tree_l1ptrel );
  reader->AddVariable("l2ptratio", &tree_l2ptratio );
  reader->AddVariable("l2ptrel", &tree_l2ptrel );
  reader->AddVariable("jet3pt", &tree_jet3pt );
  reader->AddVariable("fwd_jetpt", &tree_fwd_jetpt );

  if(ReadBDT)reader->BookMVA("BDTG method","../misc/bdt_xml/Classification_BDTG1000t2.5%n2d.weights.xml");
  //if(ReadBDT)reader->BookMVA("BDTG method","../bdts/dataset_hut_35var_loose/weights/Classification_BDTG1000t2.5%n2d.weights.xml");
  //
  TMVA::Reader *reader_hut = new TMVA::Reader( "!Color:!Silent" );
  reader_hut->AddVariable("lep1pt", &tree_lep1pt );
  reader_hut->AddVariable("lep2pt", &tree_lep2pt );
  reader_hut->AddVariable("lep1eta", &tree_lep1eta );
  reader_hut->AddVariable("lep2eta", &tree_lep2eta );
  reader_hut->AddVariable("njets", &tree_njets );
  reader_hut->AddVariable("nbtags", &tree_nbtags );
  reader_hut->AddVariable("met", &tree_met );
  reader_hut->AddVariable("ht", &tree_ht );
  reader_hut->AddVariable("mll", &tree_mll );
  reader_hut->AddVariable("nele", &tree_nele );
  reader_hut->AddVariable("jet1pt", &tree_jet1pt );
  reader_hut->AddVariable("jet2pt", &tree_jet2pt );
  reader_hut->AddVariable("btag1pt", &tree_btag1pt );
  //reader_hut->AddVariable("drl1l2", &tree_drl1l2 );
  //reader_hut->AddVariable("mindrl1j", &tree_mindrl1j );
  //reader_hut->AddVariable("mindrl2j", &tree_mindrl2j );
  //reader_hut->AddVariable("mindrl1bt", &tree_mindrl1bt );
  //reader_hut->AddVariable("mindrl2bt", &tree_mindrl2bt );
  //reader_hut->AddVariable("dphil1met", &tree_dphil1met );
  //reader_hut->AddVariable("dphil2met", &tree_dphil2met );
  reader_hut->AddVariable("mt1", &tree_mt1 );
  reader_hut->AddVariable("mt2", &tree_mt2 );
  //reader_hut->AddVariable("dphil1l2", &tree_dphil1l2 );
  // reader_hut->AddVariable("l1miniiso", &tree_l1miniiso );
  // reader_hut->AddVariable("l2miniiso", &tree_l2miniiso );
  // reader_hut->AddVariable("l1dxy", &tree_l1dxy );
  // reader_hut->AddVariable("l1dz", &tree_l1dz );
  // reader_hut->AddVariable("l2dxy", &tree_l2dxy );
  // reader_hut->AddVariable("l2dz", &tree_l2dz );
  // reader_hut->AddVariable("l1ptratio", &tree_l1ptratio );
  // reader_hut->AddVariable("l1ptrel", &tree_l1ptrel );
  // reader_hut->AddVariable("l2ptratio", &tree_l2ptratio );
  // reader_hut->AddVariable("l2ptrel", &tree_l2ptrel );
  reader_hut->AddVariable("jet3pt", &tree_jet3pt );
  reader_hut->AddVariable("fwd_jetpt", &tree_fwd_jetpt );

  //if(ReadBDT)reader_hut->BookMVA("BDTG method","../misc/bdt_xml/Classification_BDTG1000t2.5%n2d.weights_hut.xml");
  if(ReadBDT)reader_hut->BookMVA("BDTG method","../bdts/dataset_hut_NoAngleNoIso/weights/Classification_BDTG1000t2.5%n2d.weights.xml");


  //
  TMVA::Reader *reader_hut_ttbar = new TMVA::Reader( "!Color:!Silent" );
  reader_hut_ttbar->AddVariable("lep1pt", &tree_lep1pt );
  reader_hut_ttbar->AddVariable("lep2pt", &tree_lep2pt );
  reader_hut_ttbar->AddVariable("lep1eta", &tree_lep1eta );
  reader_hut_ttbar->AddVariable("lep2eta", &tree_lep2eta );
  reader_hut_ttbar->AddVariable("njets", &tree_njets );
  reader_hut_ttbar->AddVariable("nbtags", &tree_nbtags );
  reader_hut_ttbar->AddVariable("met", &tree_met );
  reader_hut_ttbar->AddVariable("ht", &tree_ht );
  reader_hut_ttbar->AddVariable("mll", &tree_mll );
  reader_hut_ttbar->AddVariable("nele", &tree_nele );
  reader_hut_ttbar->AddVariable("jet1pt", &tree_jet1pt );
  reader_hut_ttbar->AddVariable("jet2pt", &tree_jet2pt );
  reader_hut_ttbar->AddVariable("btag1pt", &tree_btag1pt );
  reader_hut_ttbar->AddVariable("mt1", &tree_mt1 );
  reader_hut_ttbar->AddVariable("mt2", &tree_mt2 );
  reader_hut_ttbar->AddVariable("jet3pt", &tree_jet3pt );
  reader_hut_ttbar->AddVariable("fwd_jetpt", &tree_fwd_jetpt );

  if(ReadBDT)reader_hut_ttbar->BookMVA("BDTG method","xmls/hut_ttbar.xml");


  //
  TMVA::Reader *reader_hut_ttv = new TMVA::Reader( "!Color:!Silent" );
  reader_hut_ttv->AddVariable("lep1pt", &tree_lep1pt );
  reader_hut_ttv->AddVariable("lep2pt", &tree_lep2pt );
  reader_hut_ttv->AddVariable("lep1eta", &tree_lep1eta );
  reader_hut_ttv->AddVariable("lep2eta", &tree_lep2eta );
  reader_hut_ttv->AddVariable("njets", &tree_njets );
  reader_hut_ttv->AddVariable("nbtags", &tree_nbtags );
  reader_hut_ttv->AddVariable("met", &tree_met );
  reader_hut_ttv->AddVariable("ht", &tree_ht );
  reader_hut_ttv->AddVariable("mll", &tree_mll );
  reader_hut_ttv->AddVariable("nele", &tree_nele );
  reader_hut_ttv->AddVariable("jet1pt", &tree_jet1pt );
  reader_hut_ttv->AddVariable("jet2pt", &tree_jet2pt );
  reader_hut_ttv->AddVariable("btag1pt", &tree_btag1pt );
  reader_hut_ttv->AddVariable("mt1", &tree_mt1 );
  reader_hut_ttv->AddVariable("mt2", &tree_mt2 );
  reader_hut_ttv->AddVariable("jet3pt", &tree_jet3pt );
  reader_hut_ttv->AddVariable("fwd_jetpt", &tree_fwd_jetpt );
  if(ReadBDT)reader_hut_ttv->BookMVA("BDTG method","xmls/hut_ttv.xml");


  
  TMVA::Reader *reader_hct = new TMVA::Reader( "!Color:!Silent" );
  reader_hct->AddVariable("lep1pt", &tree_lep1pt );
  reader_hct->AddVariable("lep2pt", &tree_lep2pt );
  reader_hct->AddVariable("lep1eta", &tree_lep1eta );
  reader_hct->AddVariable("lep2eta", &tree_lep2eta );
  reader_hct->AddVariable("njets", &tree_njets );
  reader_hct->AddVariable("nbtags", &tree_nbtags );
  reader_hct->AddVariable("met", &tree_met );
  reader_hct->AddVariable("ht", &tree_ht );
  reader_hct->AddVariable("mll", &tree_mll );
  reader_hct->AddVariable("nele", &tree_nele );
  reader_hct->AddVariable("jet1pt", &tree_jet1pt );
  reader_hct->AddVariable("jet2pt", &tree_jet2pt );
  reader_hct->AddVariable("btag1pt", &tree_btag1pt );
  // reader_hct->AddVariable("drl1l2", &tree_drl1l2 );
  // reader_hct->AddVariable("mindrl1j", &tree_mindrl1j );
  // reader_hct->AddVariable("mindrl2j", &tree_mindrl2j );
  // reader_hct->AddVariable("mindrl1bt", &tree_mindrl1bt );
  // reader_hct->AddVariable("mindrl2bt", &tree_mindrl2bt );
  // reader_hct->AddVariable("dphil1met", &tree_dphil1met );
  // reader_hct->AddVariable("dphil2met", &tree_dphil2met );
  reader_hct->AddVariable("mt1", &tree_mt1 );
  reader_hct->AddVariable("mt2", &tree_mt2 );
  // reader_hct->AddVariable("dphil1l2", &tree_dphil1l2 );
  // reader_hct->AddVariable("l1miniiso", &tree_l1miniiso );
  // reader_hct->AddVariable("l2miniiso", &tree_l2miniiso );
  // reader_hct->AddVariable("l1dxy", &tree_l1dxy );
  // reader_hct->AddVariable("l1dz", &tree_l1dz );
  // reader_hct->AddVariable("l2dxy", &tree_l2dxy );
  // reader_hct->AddVariable("l2dz", &tree_l2dz );
  // reader_hct->AddVariable("l1ptratio", &tree_l1ptratio );
  // reader_hct->AddVariable("l1ptrel", &tree_l1ptrel );
  // reader_hct->AddVariable("l2ptratio", &tree_l2ptratio );
  // reader_hct->AddVariable("l2ptrel", &tree_l2ptrel );
  reader_hct->AddVariable("jet3pt", &tree_jet3pt );
  reader_hct->AddVariable("fwd_jetpt", &tree_fwd_jetpt );

  //if(ReadBDT)reader_hct->BookMVA("BDTG method","../misc/bdt_xml/Classification_BDTG1000t2.5%n2d.weights_hct.xml");  
  if(ReadBDT)reader_hct->BookMVA("BDTG method","../bdts/dataset_hct_NoAngleNoIso/weights/Classification_BDTG1000t2.5%n2d.weights.xml");


  //
  TMVA::Reader *reader_hct_ttbar = new TMVA::Reader( "!Color:!Silent" );
  reader_hct_ttbar->AddVariable("lep1pt", &tree_lep1pt );
  reader_hct_ttbar->AddVariable("lep2pt", &tree_lep2pt );
  reader_hct_ttbar->AddVariable("lep1eta", &tree_lep1eta );
  reader_hct_ttbar->AddVariable("lep2eta", &tree_lep2eta );
  reader_hct_ttbar->AddVariable("njets", &tree_njets );
  reader_hct_ttbar->AddVariable("nbtags", &tree_nbtags );
  reader_hct_ttbar->AddVariable("met", &tree_met );
  reader_hct_ttbar->AddVariable("ht", &tree_ht );
  reader_hct_ttbar->AddVariable("mll", &tree_mll );
  reader_hct_ttbar->AddVariable("nele", &tree_nele );
  reader_hct_ttbar->AddVariable("jet1pt", &tree_jet1pt );
  reader_hct_ttbar->AddVariable("jet2pt", &tree_jet2pt );
  reader_hct_ttbar->AddVariable("btag1pt", &tree_btag1pt );
  reader_hct_ttbar->AddVariable("mt1", &tree_mt1 );
  reader_hct_ttbar->AddVariable("mt2", &tree_mt2 );
  reader_hct_ttbar->AddVariable("jet3pt", &tree_jet3pt );
  reader_hct_ttbar->AddVariable("fwd_jetpt", &tree_fwd_jetpt );

  if(ReadBDT)reader_hct_ttbar->BookMVA("BDTG method","xmls/hct_ttbar.xml");


  //
  TMVA::Reader *reader_hct_ttv = new TMVA::Reader( "!Color:!Silent" );
  reader_hct_ttv->AddVariable("lep1pt", &tree_lep1pt );
  reader_hct_ttv->AddVariable("lep2pt", &tree_lep2pt );
  reader_hct_ttv->AddVariable("lep1eta", &tree_lep1eta );
  reader_hct_ttv->AddVariable("lep2eta", &tree_lep2eta );
  reader_hct_ttv->AddVariable("njets", &tree_njets );
  reader_hct_ttv->AddVariable("nbtags", &tree_nbtags );
  reader_hct_ttv->AddVariable("met", &tree_met );
  reader_hct_ttv->AddVariable("ht", &tree_ht );
  reader_hct_ttv->AddVariable("mll", &tree_mll );
  reader_hct_ttv->AddVariable("nele", &tree_nele );
  reader_hct_ttv->AddVariable("jet1pt", &tree_jet1pt );
  reader_hct_ttv->AddVariable("jet2pt", &tree_jet2pt );
  reader_hct_ttv->AddVariable("btag1pt", &tree_btag1pt );
  reader_hct_ttv->AddVariable("mt1", &tree_mt1 );
  reader_hct_ttv->AddVariable("mt2", &tree_mt2 );
  reader_hct_ttv->AddVariable("jet3pt", &tree_jet3pt );
  reader_hct_ttv->AddVariable("fwd_jetpt", &tree_fwd_jetpt );
  if(ReadBDT)reader_hct_ttv->BookMVA("BDTG method","xmls/hct_ttv.xml");



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
      //for(unsigned int event = 0; event < 10; ++event) {
      if (STOP_REQUESTED) break;
      samesign.GetEntry(event);
      nEventsTotal++;            
      if (!quiet) bar.progress(nEventsTotal, nEventsChain);

      //Calculate weight
      //lumiAG=140;
      //lumiAG=41.4;
      //cout<<"lumiAG "<<lumiAG<<endl;
      weight = ss::is_real_data() ? 1 : ss::scale1fb()*lumiAG;
      if(proc.Contains("fcnc")) weight =  ss::scale1fb()*lumiAG*2.519*0.5/9.6; // signal samples have xsec of 9.6
      //if(proc.Contains("fcnc") and ( year == 2017 or year == 2018 )) weight =  (ss::scale1fb()*lumiAG*2.519*0.5)/(9.6*1.527); // correct for missing tau decays
      //if(proc.Contains("fcnc")) weight =  ss::scale1fb()*lumiAG*2.519*0.5;
      //if(proc.Contains("fcnc")) cout<< "fcnc weight "<<weight<<endl;
      //cout<<"weight "<<weight<<endl;      
      // Use odd events for training
      if(BDTTraining and event %2 == 0) continue ;
      //Use even events for Application
      if(BDTApplication and event %2 == 1) continue ;
      if(BDTTraining or BDTApplication) weight = weight*2;
      
      // Simple cuts first to speed things up
      lep1ccpt = ss::lep1_coneCorrPt();
      lep2ccpt = ss::lep2_coneCorrPt();
      if (not noLeptonPtCut) {
	if (lep1ccpt < 25) continue;
	if (lep2ccpt < 20) continue;
      }
      lep1id = ss::lep1_id();
      lep2id = ss::lep2_id();
      //bool pass_trig = (doSS) ? ss::fired_trigger_ss() : ss::fired_trigger();	    

      bool pass_trig =  ss::fired_trigger_ss();
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
      //nleps = ss::nleps();

      njets40 = 0; float ht40 = 0; 
      if(njets>0)for(int i =0; i< njets;i++){
	  if(ss::jets()[i].pt()>40){
	    njets40++;
	    ht40 = ht40+ss::jets()[i].pt();
	  }
	}
      nbtags25 = 0;      nbtagsM = 0;       nbtags25M = 0;
      if(nbtags>0)for(int i=0;i<nbtags;i++){
	  if(ss::btags()[i].pt()>25) nbtags25++;
	  nbtagsM++;
	  if(ss::btags()[i].pt()>25) nbtags25M++;	  	  
	}
      
      
      nleps = (lep3good) ? ((ss::lep4_passes_id() and (ss::lep4_p4().pt() > (abs(ss::lep4_id())==11 ? 15 : 10))) ? 4 : 3) : 2;      
      vector<lepton> lep;
      lepton temp;
      temp.v = ss::lep1_p4(); temp.id = ss::lep1_id(); temp.isgood = lep1good;
      lep.push_back(temp);
      temp.v = ss::lep2_p4(); temp.id = ss::lep2_id(); temp.isgood = lep2good;
      lep.push_back(temp);
      temp.v = ss::lep3_p4(); temp.id = ss::lep3_id(); temp.isgood = lep3good;
      if(lep3good)lep.push_back(temp);
      //if(MOSSF(lep)<999999)cout<<"mossf "<<MOSSF(lep)<<endl; 
      
      bool isClass6 = ss::hyp_class() == 6;
      float mtnonz = ss::mtmin();
      SR = signal_region(ss::njets(), ss::nbtags(), ss::met(), ss::ht(), ss::mtmin(), lep1id, lep2id, lep1ccpt, lep2ccpt, lep3ccpt, nleps);
      //SR = signal_region(njets40, nbtags25M, ss::met(), ht40 , ss::mtmin(), lep1id, lep2id, lep1ccpt, lep2ccpt, lep3ccpt, nleps); 
      ht = ss::ht();
      nisrmatch = ss::nisrMatch();      	    
      /* hyp_class
       * 1: SS, loose-loose
       * 2: SS, tight-loose (or loose-tight)
       * 3: SS, tight-tight
       * 4: OS, tight-tight
       * 5: SS, inSituFR
       * 6: SS, tight-tight and fails Z-veto (lies! hyp_class==6 != lep1good and lep2good)
       */
      int hyp_class = ss::hyp_class();      
      if (!ss::is_real_data()) {
	weight *= getTruePUw(year, ss::trueNumInt()[0]);
	if (lep1good) weight *= leptonScaleFactor(year, lep1id, lep1ccpt, lep1eta, ht);
	if (lep2good) weight *= leptonScaleFactor(year, lep2id, lep2ccpt, lep2eta, ht);
	if (lep3good) weight *= leptonScaleFactor(year, lep3id, lep3ccpt, lep3eta, ht);

	if (!lep3good) {
	  weight *= triggerScaleFactor(year, lep1id, lep2id, lep1pt, lep2pt, lep1eta, lep2eta, ht, analysis, 0);
	}
	weight *= ss::weight_btagsf();
	if (proc.Contains("ttw"))  weight *= isrWeight(year, nisrmatch, 1);
	if (proc.Contains("ttz"))  weight *= isrWeight(year, nisrmatch, 2);
	if (proc.Contains("tt_"))  weight *= isrWeight(year, nisrmatch, 10);

	if (year == 2016) weight *= ss::prefire2016_sf();
	if (year == 2017) weight *= ss::prefire2017_sf();		
	if (useTTBB and (ss::extragenb() >= 2)) weight *= 1.7;	
	weight *=ss::decayWSF();	
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
	  // subtract prompt MC
	  if (isData and !ss::is_real_data()) weight *= -1.;
	  hyp_class = 3; // we've faked a SS Tight-Tight with a SS LL or SS TL
	  // Basically just update this so it gets put in the SR	  	  
	} else {
	  continue; // Not a fakeing hyp_class
	}
      }      

      //
      // bool class6Fake = false;     
      // if (doFakes) {
      // 	if (hyp_class == 1 or hyp_class == 2) {
      // 	  bool foundGoodLoose = false;
      // 	  if (ss::lep1_passes_id()==0 ) {
      // 	    float fr = fakeRate(year, lep1id, lep1ccpt, lep1eta, ht, analysis, new2016FRBins, !minPtFake18);
      // 	    weight *= fr/(1.-fr);
      // 	    foundGoodLoose = true;
      // 	  }
      // 	  if (ss::lep2_passes_id()==0 ) {
      // 	    float fr = fakeRate(year, lep2id, lep2ccpt, lep2eta, ht, analysis, new2016FRBins, !minPtFake18);
      // 	    weight *= fr/(1.-fr);
      // 	    foundGoodLoose = true;
      // 	  }
      // 	  if (!foundGoodLoose)
      // 	    continue;
      // 	  // subtract double FO (why is this?)                                                                                                                                      
      // 	  if (hyp_class == 1) weight *= -1.;
      // 	  hyp_class = 3; // we've faked a SS Tight-Tight with a SS LL or SS TL                                                                                                      
      // 	  // Basically just update this so it gets put in the SR                                                                                                     
      // 	} 

      // 	else {
      // 	  continue; // Not a fakeing hyp_class                                                                                                                                      
      // 	}
      // }

      //flips

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

      // if all 3 charges are the same, throw the event away
      if (nleps > 2 and ((lep1id>0 and lep2id>0 and lep3id>0) or
			 (lep1id<0 and lep2id<0 and lep3id<0))) continue;
      

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
	    
      // jet pt 
      // ptj1 = (njets >= 1) ? ss::jets()[0].pt()*ss::jets_undoJEC()[0]*ss::jets_JEC()[0] : 0;
      // ptj2 = (njets >= 2) ? ss::jets()[1].pt()*ss::jets_undoJEC()[1]*ss::jets_JEC()[1] : 0;
      // ptj3 = (njets >= 3) ? ss::jets()[2].pt()*ss::jets_undoJEC()[2]*ss::jets_JEC()[2] : 0;
      // ptj4 = (njets >= 4) ? ss::jets()[3].pt()*ss::jets_undoJEC()[3]*ss::jets_JEC()[3] : 0;

      // ptbt1 = (nbtags >= 1) ? ss::btags()[0].pt()*ss::btags_undoJEC()[0]*ss::btags_JEC()[0] : 0;
      // ptbt2 = (nbtags >= 2) ? ss::btags()[1].pt()*ss::btags_undoJEC()[1]*ss::btags_JEC()[1] : 0;
      // ptbt3 = (nbtags >= 3) ? ss::btags()[2].pt()*ss::btags_undoJEC()[2]*ss::btags_JEC()[2] : 0;
      // ptbt4 = (nbtags >= 4) ? ss::btags()[3].pt()*ss::btags_undoJEC()[3]*ss::btags_JEC()[3] : 0;

      
      ptj1 = (njets >= 1) ? ss::jets()[0].pt() : 0;
      ptj2 = (njets >= 2) ? ss::jets()[1].pt() : 0;
      ptj3 = (njets >= 3) ? ss::jets()[2].pt() : 0;
      ptj4 = (njets >= 4) ? ss::jets()[3].pt() : 0;

      ptbt1 = (nbtags >= 1) ? ss::btags()[0].pt() : 0;
      ptbt2 = (nbtags >= 2) ? ss::btags()[1].pt() : 0;
      ptbt3 = (nbtags >= 3) ? ss::btags()[2].pt() : 0;
      ptbt4 = (nbtags >= 4) ? ss::btags()[3].pt() : 0;

      
      // float ht_new  = 0;
      // for(int i =0; i<njets;i++){
      // 	ht_new  = ht_new + ss::jets()[i].pt();
      // }

      // float ht_new2  = 0;
      // for(int i =0; i<njets;i++){
      // 	ht_new2  = ht_new2 + ss::jets()[i].pt()*ss::jets_undoJEC()[i]*ss::jets_JEC()[i];
      // }
      
      // if(njets >= 1)cout<<" ptj1 "<<ptj1<< " "<<ss::jets()[0].pt()<<endl;
      // cout<<" ht "<<ht<<" "<<ht_new<<" "<<ht_new2<<endl;
     
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
			   do_fill(h_mll, m12);			   
			   do_fill(h_mllzoom, m12);
			   do_fill(h_mossf,MOSSF(lep));			   
			   //if (nleps > 2) do_fill(h_zmll, mllos);
			   if (nleps > 2) do_fill(h_m3l, m3l);
			   do_fill(h_nleps, nleps);
			   do_fill(h_mtmin, ss::mtmin());
			   do_fill(h_dphil1met, calcDeltaPhi(lep1phi,metphi));
			   do_fill(h_dphil2met, calcDeltaPhi(lep2phi,metphi));
			   if(njets>0)do_fill(h_dphimetj1, calcDeltaPhi(ss::jets()[0].phi(),metphi));
			   do_fill(h_dphil1l2, calcDeltaPhi(lep1phi,lep2phi));	
			   do_fill(h_njets, njets);
			   do_fill(h_njets40, njets40);
			   do_fill(h_nisrjets, nisrjets);
			   do_fill(h_nisrmatch, nisrmatch);	        
			   do_fill(h_ntb40, ntb40);
			   do_fill(h_nbtags, nbtags);
			   do_fill(h_nbtagszoom, nbtags);
			   do_fill(h_nbtags25, nbtags25);
			   do_fill(h_nbtagsM, nbtagsM);
			   do_fill(h_nbtags25M, nbtags25M);	
			   //int nbnj = 5*min(nbtags,3)+(max(min(njets,6),2)-2);
			   //do_fill(h_nbnj, nbnj);
			   do_fill(h_pt1, lep1ccpt);
			   do_fill(h_pt2, lep2ccpt);
			   if (nleps > 2) do_fill(h_pt3, lep3pt);
			   do_fill(h_eta1,   ss::lep1_p4().eta());
			   do_fill(h_eta2,   ss::lep2_p4().eta());

			   int looseleg = -1;
			   if (hyp_class == 2) {
			     looseleg = (lep1good ? 2 : 1);
			   }
			   // if (looseleg == 1) do_fill(abs(lep1id) == 11 ? h_etaelnt     : h_etamlnt,     lep1eta);
			   // if (looseleg == 2) do_fill(abs(lep2id) == 11 ? h_etaelnt     : h_etamlnt,     lep2eta);
			   do_fill(h_ptrel1, lep1ptrel);
			   do_fill(h_ptrel2, lep2ptrel);
			   //if (looseleg > 0) do_fill(h_ptrellnt, looseleg == 1 ? lep1ptrel   : lep2ptrel);
			   do_fill(abs(lep1id) == 11 ? h_ptrele   : h_ptrelm,   lep1ptrel);
			   do_fill(abs(lep2id) == 11 ? h_ptrele   : h_ptrelm,   lep2ptrel);
			   // do_fill2D(abs(lep1id) == 11 ? h_ptabsetae     : h_ptabsetam,     lep1pt, fabs(lep1eta));
			   // do_fill2D(abs(lep2id) == 11 ? h_ptabsetae     : h_ptabsetam,     lep2pt, fabs(lep2eta));
			   do_fill(h_ptratio1, lep1ptratio);
			   do_fill(h_ptratio2, lep2ptratio);
			   do_fill(h_miniiso1, lep1miniiso);
			   do_fill(h_miniiso2, lep2miniiso);	
			   int type = ss::hyp_type();
			   do_fill(h_type,   type>1 ? type-1 : type);
			   //do_fill(h_q1,   lep1id>0 ? -1 : 1);
			   if (nleps > 2) {
			     if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 39) do_fill(h_type3l, 0); // mu mu mu
			     if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 37) do_fill(h_type3l, 1); // mu mu e
			     if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 35) do_fill(h_type3l, 2); // mu e e
			     if (abs(lep1id) + abs(lep2id) + abs(lep3id) == 33) do_fill(h_type3l, 3); // e e e
			   }
			   do_fill(h_nvtx,   ss::nGoodVertices());
			   do_fill(h_sr,   SR);
			   if(ReadBDT) do_fill(h_bdt,  reader->EvaluateMVA("BDTG method"));
			   if(ReadBDT) do_fill(h_bdt_hut,  reader_hut->EvaluateMVA("BDTG method"));
			   if(ReadBDT) do_fill(h_bdt_hct,  reader_hct->EvaluateMVA("BDTG method"));

			   if(ReadBDT) do_fill(h_bdt_hut_ttbar,  reader_hut_ttbar->EvaluateMVA("BDTG method"));
			   if(ReadBDT) do_fill(h_bdt_hut_ttv,  reader_hut_ttv->EvaluateMVA("BDTG method"));
			   if(ReadBDT) do_fill2D(h_bdt2D_hut,  reader_hut_ttbar->EvaluateMVA("BDTG method"), reader_hut_ttv->EvaluateMVA("BDTG method"));
			   			  			   
			   if(ReadBDT) do_fill(h_bdt_hct_ttbar,  reader_hct_ttbar->EvaluateMVA("BDTG method"));
			   if(ReadBDT) do_fill(h_bdt_hct_ttv,  reader_hct_ttv->EvaluateMVA("BDTG method"));
			   if(ReadBDT) do_fill2D(h_bdt2D_hct,  reader_hct_ttbar->EvaluateMVA("BDTG method"), reader_hct_ttv->EvaluateMVA("BDTG method"));
			   
			   
			   do_fill(h_drl1l2, calcDeltaR(ss::lep1_p4().eta(),ss::lep1_p4().phi(),ss::lep2_p4().eta(),ss::lep2_p4().phi()));
			   do_fill(h_mindrl1j, minDR(ss::lep1_p4(),ss::jets()));
			   do_fill(h_mindrl2j, minDR(ss::lep2_p4(),ss::jets()));
			   do_fill(h_mindrl1bt, minDR(ss::lep1_p4(),ss::btags()));
			   do_fill(h_mindrl2bt, minDR(ss::lep2_p4(),ss::btags()));
			   do_fill(h_mt1,calcMT(ss::lep1_p4().pt(),ss::lep1_p4().phi(), met, metphi));
			   do_fill(h_mt2,calcMT(ss::lep2_p4().pt(),ss::lep2_p4().phi(), met, metphi));
			   do_fill(h_l1dxy, ss::lep1_dxyPV());
			   do_fill(h_l2dxy, ss::lep2_dxyPV());
			   do_fill(h_l1dz,ss::lep1_dZ());
			   do_fill(h_l2dz,ss::lep2_dZ());
			   do_fill(h_ptj1,   ptj1);
			   do_fill(h_ptj2,   ptj2);
			   do_fill(h_ptj3,   ptj3);
			   do_fill(h_ptj4,   ptj4);
			   do_fill(h_ptbt1,   ptbt1);
			   do_fill(h_ptbt2,   ptbt2);
			   do_fill(h_ptbt3,   ptbt3);
			   do_fill(h_ptbt4,   ptbt4);
			   do_fill(h_fwd_jetpt, PtMaxEta(ss::jets()));
			   do_fill(h_nele,  (abs(lep1id)==11 and abs(lep2id)==11) ? 2 : (abs(lep1id)==11 or abs(lep2id)==11) ? 1 : 0  );
			 };

      bool BR_lite = ht > 300  and njets >= 2 and  nbtags >= 2 and met >= 50;
      bool BR = BR_lite and hyp_class == 3;
      if (BR) fill_region("br", weight);
      if (hyp_class == 2) fill_region("susytl", weight);
      if (hyp_class == 2 and lep1ccpt > 25 and lep2ccpt > 25) fill_region("hhtl", weight);            	    
      if (hyp_class == 3 and nleps == 2 and njets >= 2 and met > 50 and lep1ccpt > 25 and lep2ccpt > 25 ) {
	//if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and met > 50 and lep1ccpt > 25 and lep2ccpt > 25 ) {	      	
	fill_region("sshh", weight);
      }
      bool OnZ = fabs(m12-91.2) <15.;            
      bool LowMetOnZor0b = (met<50&&OnZ)||nbtags==0;
      bool MossfOnZ = fabs(MOSSF(lep)-91.2)<15.;
      bool mllowmetonz2b = met<50. and nleps >2 and hyp_class == 6 and nbtags >= 2 and  lep1good and lep2good and lep3good and MossfOnZ;         
      //if  ((hyp_class == 1 or hyp_class == 2 or hyp_class == 3 ) and nleps == 2 and njets >= 2 and nbtags > 0 and lep1ccpt > 25 and lep2ccpt > 20 ){ 
      if  ( hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags > 0 and lep1ccpt > 25 and lep2ccpt > 20 ){ 
	fill_region("ssbr", weight); 
	//fill the tree variables
	tree_lep1pt = lep1pt;
	tree_lep2pt = lep2pt;
	tree_lep1eta = lep1eta;
	tree_lep2eta = lep2eta;
	tree_njets = njets;
	tree_nbtags = nbtags;
	tree_met = met;
	tree_ht = ht;
	tree_mll = m12;
	tree_nele = (abs(lep1id)==11 and abs(lep2id)==11) ? 2 : (abs(lep1id)==11 or abs(lep2id)==11) ? 1 : 0 ;
	tree_jet1pt = ptj1;
	tree_jet2pt = ptj2;
	tree_btag1pt = ptbt1;
	tree_drl1l2 =  calcDeltaR(ss::lep1_p4().eta(),ss::lep1_p4().phi(),ss::lep2_p4().eta(),ss::lep2_p4().phi());
	tree_mindrl1j = minDR(ss::lep1_p4(),ss::jets());
	tree_mindrl2j = minDR(ss::lep2_p4(),ss::jets());
	tree_mindrl1bt = minDR(ss::lep1_p4(),ss::btags());
	tree_mindrl2bt = minDR(ss::lep2_p4(),ss::btags());
	tree_dphil1met = calcDeltaPhi(lep1phi,metphi);
	tree_dphil2met = calcDeltaPhi(lep2phi,metphi);
	tree_mt1 = calcMT(ss::lep1_p4().pt(),ss::lep1_p4().phi(), met, metphi);
	tree_mt2 = calcMT(ss::lep2_p4().pt(),ss::lep2_p4().phi(), met, metphi);
	tree_dphil1l2 = calcDeltaPhi(lep1phi,lep2phi);
	tree_l1miniiso = ss::lep1_miniIso();
	tree_l2miniiso = ss::lep2_miniIso();
	tree_l1dxy = ss::lep1_dxyPV();
	tree_l1dz = ss::lep1_dZ();
	tree_l2dxy = ss::lep2_dxyPV();
	tree_l2dz = ss::lep2_dZ();
	tree_l1ptratio = lep1ptratio;
	tree_l1ptrel = lep1ptrel;
	tree_l2ptratio = lep2ptratio;
	tree_l2ptrel = lep2ptrel;
	tree_jet3pt = ((int)ss::jets().size() >= 3) ? ss::jets()[2].pt() : 0;
        tree_fwd_jetpt = PtMaxEta(ss::jets());	
	tree_weight = weight;	
	out_tree->Fill();	
	//cout<<" bdt value"<<reader_hut->EvaluateMVA("BDTG method")<<endl;
      }
      if (hyp_class == 3 and nleps == 2 and met < 50 and lep1ccpt > 25 and lep2ccpt > 20 and abs(lep1id) == 11 and abs(lep1id) == 11 and abs(m12-91.2) <15. ) fill_region("ssonz", weight);
      //if (hyp_class == 3 and nleps == 2 and njets >= 2 and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ssbr", weight);      
      if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags == 0  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss0b2j", weight); 
      if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2j", weight); 
      if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2j", weight);
      
      // if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtagsM == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2jbtagM", weight); 
      // if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtagsM  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2jbtagM", weight);
      
      // if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags25 == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2jbtag25", weight); 
      // if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags25  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2jbtag25", weight);
      
      // if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags25M == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2jbtag25M", weight); 
      // if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags25M  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2jbtag25M", weight);

      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtags == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2jjet40", weight); 
      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtags  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2jjet40", weight);
      
      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtagsM == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2jjet40btagM", weight); 
      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtagsM  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2jjet40btagM", weight);
      
      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtags25 == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2jjet40btag25", weight); 
      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtags25  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2jjet40btag25", weight);
      
      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtags25M == 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss1b2jjet40btag25M", weight); 
      // if (hyp_class == 3 and nleps == 2 and njets40 >= 2 and nbtags25M  > 1  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("ss2b2jjet40btag25M", weight);
           
      if (hyp_class == 3 and LowMetOnZor0b and lep1ccpt > 25 and lep2ccpt > 20 )  fill_region("lowmetonzor0b", weight);	      
      if (hyp_class == 3 and nleps == 2 and njets >= 2 and lep1ccpt > 25 and lep2ccpt > 20 and !LowMetOnZor0b ) fill_region("ssbr2", weight); 
      if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags == 1   and lep1ccpt > 25 and lep2ccpt > 20 and !LowMetOnZor0b ) fill_region("ss1b2j2", weight); 
      if (hyp_class == 3 and nleps == 2 and njets >= 2 and nbtags  > 1   and lep1ccpt > 25 and lep2ccpt > 20 and !LowMetOnZor0b ) fill_region("ss2b2j2", weight);      
      
      if (hyp_class == 3 and nleps >  2 and njets >= 1 and nbtags >  0 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 ) fill_region("mlbr", weight);  
      if (hyp_class == 3 and nleps >  2 and njets >= 1 and nbtags == 1 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 ) fill_region("ml1b1j", weight);  
      if (hyp_class == 3 and nleps >  2 and njets >= 1 and nbtags >  1 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 ) fill_region("ml2b2j", weight);
      
      if (mllowmetonz2b  and njets >= 1 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 ) fill_region("mllowmetonz2b", weight);  

      if ((hyp_class == 3 or hyp_class == 6) and nleps >  2 and njets >= 1 and nbtags >  0 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 and !mllowmetonz2b ) fill_region("mlbrinc", weight);  
      if ((hyp_class == 3 or hyp_class == 6) and nleps >  2 and njets >= 1 and nbtags == 1 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 and !mllowmetonz2b ) fill_region("ml1b1jinc", weight);  
      if ((hyp_class == 3 or hyp_class == 6) and nleps >  2 and njets >= 1 and nbtags >  1 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 and !mllowmetonz2b ) fill_region("ml2b2jinc", weight);
      
      if (hyp_class == 6 and nleps >  2 and njets >= 1 and nbtags >  0 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 and !mllowmetonz2b ) fill_region("mlbronz", weight);  
      if (hyp_class == 6 and nleps >  2 and njets >= 1 and nbtags == 1 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 and !mllowmetonz2b ) fill_region("ml1b1jonz", weight);  
      if (hyp_class == 6 and nleps >  2 and njets >= 1 and nbtags >  1 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 and !mllowmetonz2b ) fill_region("ml2b2jonz", weight);
      
      //if (hyp_class == 3 and nleps >  2 and njets >= 2 and nbtags == 2 and lep1ccpt > 25 and lep2ccpt > 20 and lep3ccpt > 20 and MossfOnZ ) fill_region("ttzbr", weight);
      if (hyp_class == 4 and nleps == 2 and njets >= 2 and nbtags > 0  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("osbr", weight);
      if (hyp_class == 2 and nleps == 2 and njets >= 2 and nbtags > 0  and lep1ccpt > 25 and lep2ccpt > 20 ) fill_region("tl", weight); 
      
    }//event loop
    //delete tree;
    delete file;
  }//file loop

  f1->cd();
  for (HistCol* coll : registry) coll->Write();
  for (HistCol2D* coll : registry2D) coll->Write();
  out_tree->Write();
  f1->Close();
  if (!quiet) cout << "\n Done!" << endl;
  return 0;
}

