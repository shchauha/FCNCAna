#ifndef SSBABYMAKER_H
#define SSBABYMAKER_H

#include <vector>
#include <functional>
#include <algorithm>
#include <chrono>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH2.h"
#include "Math/VectorUtil.h"
#include "CORE/CMS3.h"
#include "CORE/Config.h"
#include "CORE/Tools/utils.h"
#include "CORE/Tools/JetCorrector.h"
#include "CORE/Tools/MT2/MT2.h"
#include "CORE/SSSelections.h"
#include "CORE/ElectronSelections.h"
#include "CORE/IsolationTools.h"
#include "CORE/JetSelections.h"
#include "CORE/MetSelections.h"
#include "CORE/MCSelections.h"
#include "CORE/Tools/goodrun.h"
#include "CORE/Tools/JetResolution/JetResolution.h"
#include "CORE/Tools/jetcorr/JetCorrectionUncertainty.h"
#include "CORE/Tools/jetcorr/SimpleJetCorrectionUncertainty.h"
#include "CORE/Tools/datasetinfo/getDatasetInfo.h"
#include "TROOT.h"
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "fastsim/xsec.h"


#ifdef __MAKECINT__
#pragma link C++ class ROOT::Math::PxPyPzE4D<float>+;
#pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >+;
#pragma link C++ typedef ROOT::Math::XYZTVectorF;
#endif
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"

#ifdef __MAKECINT__
#pragma link C++ class ROOT::Math::PxPyPzE4D<float>+;
#pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >+;
#pragma link C++ typedef ROOT::Math::XYZTVectorF;
#endif
#include "Math/Vector4D.h"
#include "Math/LorentzVector.h"

#ifdef __MAKECINT__
#pragma link C++ class ROOT::Math::PxPyPzE4D<float>+;
#pragma link C++ class ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >+;
#pragma link C++ typedef ROOT::Math::XYZTVectorF;
#endif

#define PI 3.14159265


int nPoints_higgs(int sample, int mH);
float xsec_higgs(int sample, int mH);

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef vector<pair<const LorentzVector *, double> > jets_with_corr_t;

class BTagCalibration;
class BTagCalibrationReader;

using namespace std;
using namespace tas;

struct csErr_t { float cs_scale_no = 0; float cs_scale_up = 0; float cs_scale_dn = 0; float cs_pdf[102] = {0}; int SR = -1; bool isGood = 0; };

//Classes
class babyMaker {

  public:
    babyMaker(bool debug = 0) {
      path = ".";
      verbose = debug;
      evt_cut = 0;
    }
    void MakeBabyNtuple(const char* output_name, bool isFastsim, int iSignal);
    void InitBabyNtuple();
    void CloseBabyNtuple () {
        BabyFile->cd();
        BabyTree->Write();
        BabyFile->Close();
    }
    csErr_t ProcessBaby(string filename_in, FactorizedJetCorrector* jetCorr, JetCorrectionUncertainty *jetUnc, bool isFastsim = 0, int iSignal = 0);

    bool ignore_scale1fb = false;
    bool ignore_os = false;
    bool read_psweights = false;
    bool has_lhe_branches = false;
    bool verbose = false;
    unsigned int evt_cut = 0;

    TFile* BabyFile;
    TTree* BabyTree;

  protected:
    //

  private:

    //Switches
    TString path;


    // for btag SFs
    float getBtagEffFromFile(float pt, float eta, int mcFlavour);

    BTagCalibration* calib1;
    BTagCalibration* calib2;
    BTagCalibration* calib3;
    BTagCalibration* calib4;
    BTagCalibration* calib_fs;

    BTagCalibrationReader* btcr1;
    BTagCalibrationReader* btcr2;
    BTagCalibrationReader* btcr3;
    BTagCalibrationReader* btcr4;
    BTagCalibrationReader* btcr_fs;
    BTagCalibrationReader* btcr_test;

    TH2D* h_btag_eff_b;
    TH2D* h_btag_eff_c;
    TH2D* h_btag_eff_udsg;

    // JER
  JetResolution res;  

    // DatasetInfo
  DatasetInfoFromFile df;


    //MET
    float met;
    float metPhi;
    float rawmet;
    float rawmetPhi;
    float calomet;
    float calometPhi;

    //Meta Variables
    ULong64_t event;
    int lumi;
    int run;
    bool is_real_data;
    bool is_fastsim;
    bool is_miniaodv1_80X;
    float scale1fb;
    float genps_weight;
    float xsec;
    int neventstotal;
    float xsec_error;
    float xsec_ps;
    vector <TString> sparmNames;
    vector <float> sparms;
    float kfactor;
    string filename;
    // vector <float> genweights;
    vector <float> pdfweights;
    vector <string> genweightsID;
    bool passedFilterList;
    float prefire_sf;
    float prefire_sfup;
    float prefire_sfdown;
    float prefire2016_sf;
    float prefire2016_sfup;
    float prefire2016_sfdown;
    float prefire2017_sf;
    float prefire2017_sfup;
    float prefire2017_sfdown;
    float prefire2017ele_sf;
    float prefire2017ele_sfup;
    float prefire2017ele_sfdown;

    //Pileup
    vector <float> trueNumInt;
    vector <int> nPUvertices;
    int nGoodVertices;

    //Filters
    bool filt_csc;
    bool filt_hbhe;
    bool filt_hcallaser;
    bool filt_ecaltp;
    bool filt_trkfail;

    //Gen MET
    float gen_met;
    float gen_met_phi;

    bool skim;
    bool skim_nomet;
    bool br;
    int sr;
    int nleps;

    int ss_nleps;
    bool ss_br;
    int ss_region;
    int ss_sr;

    int extragenb;
    int extragenc;
    int ngenjets;
    float genht;
    int ngenjets30;
    float genht30;
    bool passfilter;
    vector <int> extragenbmoms;
    vector <int> extragencmoms;

    bool hasgammatoll;
    bool gammatollmomemu;
    float gammatolldr;
    LorentzVector gammatoll1;
    LorentzVector gammatoll2;

    // for TTTT
    vector <int> bjet_type;
    vector <int> jet_type;
    int ndrlt0p4;
    bool gengood;
    int nleptonic;
    int ntau;
    int nleptonicW;
    float mfourtop;
    int nhadronicW;
    int nW;
    float leptonicWSF;
    float hadronicWSF;
    float decayWSF;
    int njincone1;
    int njincone2;

    //Jets
    int njets_raw;
    float ht_raw;
    int njets;
    int njetsAG;
    int nbtagsAG;
    float ht;
    vector <LorentzVector> jets;
    vector <LorentzVector> jets_jec_up;
    vector <LorentzVector> jets_jec_dn;
    vector <LorentzVector> bjets_jec_up;
    vector <LorentzVector> bjets_jec_dn;
    vector <float> jets_disc_up;
    vector <float> jets_disc_dn;
    vector <int> jets_flavor;
    vector <float> jets_disc;
    vector <float> jets_bsf;
    vector <float> jets_disc_mva;
    vector <float> jets_disc_ivf;
    vector <float> jets_JEC;
    vector <float> jets_unc;
    vector <float> jets_undoJEC;
    vector <LorentzVector> mostJets; //this is all jets with pT > 5, |eta| < 2.4
    vector <LorentzVector> mostJets_rawp4; // this is all jets with pT > 5, |eta| < 2.4
    vector <float> mostJets_disc;          // this is all jets with pT > 5, |eta| < 2.4
    vector <float> mostJets_unc;           // this is all jets with pT > 5, |eta| < 2.4
    vector <float> mostJets_JEC;           // this is all jets with pT > 5, |eta| < 2.4
    vector <float> mostJets_undoJEC;       // this is all jets with pT > 5, |eta| < 2.4
    vector <bool> mostJets_passCleaning;
    vector <bool> mostJets_passID;
    // vector <float> jet_pt;

    //Jets unc stuff
    int njets_unc_up;
    int njets_unc_dn;
    float ht_unc_up;
    float ht_unc_dn;
    int nbtags_unc_up;
    int nbtags_unc_dn;
    float met_unc_up;
    float met_unc_dn;
    float metPhi_unc_up;
    float metPhi_unc_dn;

    //Jets unc stuff - JER
    int njets_JER;
    int njets_JER_up;
    int njets_JER_dn;
    int nbtags_JER;
    int nbtags_JER_up;
    int nbtags_JER_dn;
    float ht_JER;
    float ht_JER_up;
    float ht_JER_dn;
    float met_JER;
    float met_JER_up;
    float met_JER_dn;
    float metPhi_JER;
    float metPhi_JER_up;
    float metPhi_JER_dn;

    int hyp_class;

    //Gen flags
    bool lep1_isPrompt;
    bool lep1_isDirectPrompt;
    bool lep1_isStat3;
    bool lep2_isPrompt;
    bool lep2_isDirectPrompt;
    bool lep2_isStat3;
    bool lep2_genps_isHardProcess;
    bool lep2_genps_fromHardProcessFinalState;
    bool lep2_genps_fromHardProcessDecayed;
    bool lep2_genps_isDirectHardProcessTauDecayProductFinalState;
    bool lep2_genps_fromHardProcessBeforeFSR;
    bool lep2_genps_isLastCopy;
    bool lep2_genps_isLastCopyBeforeFSR;
    bool lep1_genps_isHardProcess;
    bool lep1_genps_fromHardProcessFinalState;
    bool lep1_genps_fromHardProcessDecayed;
    bool lep1_genps_isDirectHardProcessTauDecayProductFinalState;
    bool lep1_genps_fromHardProcessBeforeFSR;
    bool lep1_genps_isLastCopy;
    bool lep1_genps_isLastCopyBeforeFSR;
    int lep1_mc3idx = -1;
    int lep2_mc3idx = -1;

    //Leptons
    LorentzVector lep1_p4;
    LorentzVector lep2_p4;
    LorentzVector dilep_p4;
    float dilep_mass;
    float mass13;
    float mass23;
    float mass123;
    int lep1_id;
    int lep2_id;
    int lep1_idx;
    int lep2_idx;
    int hyp_type;
    float lep1_coneCorrPt;
    float lep2_coneCorrPt;
    float lep3_coneCorrPt;
    float lep4_coneCorrPt;
    float lep1_pt;
    float lep2_pt;
    float lep3_pt;
    float lep1_eta;
    float lep2_eta;
    float lep3_eta;
    float lep1_phi;
    float lep2_phi;
    float lep3_phi;

    //Lepton Mother (old isFromW function)
    int lep1_motherID;
    int lep2_motherID;

    int lep1_mc_id;
    int lep2_mc_id;
    int lep1_mc_motherid;
    int lep2_mc_motherid;

    //Electron-specific
    bool lep1_el_conv_vtx_flag;
    int lep1_el_exp_innerlayers;
    bool lep1_nPixelMiss;
    bool lep1_tightCharge;
    bool lep1_mu_ptErr;
    bool lep1_el_threeChargeAgree;
    float lep1_el_etaSC;
    bool lep1_isTrigSafeNoIsov1;
    bool lep1_isTrigSafev1;

    bool lep2_el_conv_vtx_flag;
    int lep2_el_exp_innerlayers;
    bool lep2_nPixelMiss;
    bool lep2_tightCharge;
    bool lep2_mu_ptErr;
    bool lep2_el_threeChargeAgree;
    float lep2_el_etaSC;
    bool lep2_isTrigSafeNoIsov1;
    bool lep2_isTrigSafev1;

    // Muon-specific
    bool lep1_tightMuonPOG;
    bool lep1_mediumMuonPOG;
    bool lep1_looseMuonPOG;
    bool lep2_tightMuonPOG;
    bool lep2_mediumMuonPOG;
    bool lep2_looseMuonPOG;
    bool lep3_mediumMuonPOG;

    // FG: add VVV lepton info
    bool lep1_el_vvv_ss_fo_id;
    bool lep1_el_vvv_3l_fo_id;
    bool lep1_el_vvv_ss_id;
    bool lep1_el_vvv_3l_id;
    bool lep1_mu_vvv_fo_id;
    bool lep1_mu_vvv_id;
    float lep1_vvv_iso;
    bool lep2_el_vvv_ss_fo_id;
    bool lep2_el_vvv_3l_fo_id;
    bool lep2_el_vvv_ss_id;
    bool lep2_el_vvv_3l_id;
    bool lep2_mu_vvv_fo_id;
    bool lep2_mu_vvv_id;
    float lep2_vvv_iso;

    //b-tags
    int nbtags_raw;
    vector <LorentzVector> btags;
    vector <int> btags_flavor;
    vector <float> btags_disc;
    vector <float> btags_cdisc;
    vector <float> btags_disc_mva;
    vector <float> btags_disc_ivf;
    vector <float> btags_JEC;
    vector <float> btags_unc;
    vector <float> btags_eff;
    vector <float> btags_effpt;
    vector <float> btags_sf;
    vector <float> btags_undoJEC;
    int nbtags;

    //Scale factors (from 8 TeV, outdated)
    float sf_dilepTrig_hpt;
    float sf_dilepTrig_lpt;
    float sf_dilepTrig_vlpt;
    float sf_dilep_eff;

    //mT for both leptons, mt2
    float mt;
    float mt_l2;
    float mt2;
    float mtmin;

    //Look at gen-level particles and choose favorite hypothesis
    //Limited usefulness as-is; probably not smart enough to reject non-prompt leptons
    int lep1_id_gen;
    int lep2_id_gen;
    LorentzVector lep1_p4_gen;
    LorentzVector lep2_p4_gen;

    //Third lepton -- quality is 2 for good, 1 for fakeable, 0 for veto
    int lep3_id;
    int lep3_idx;
    int lep3_mcid;
    int lep3_mc_motherid;
    int lep3_motherID;
    int lep3_mc3idx;
    int lep3_quality;
    LorentzVector lep3_p4;

    ////Fourth lepton -- quality is 2 for good, 1 for fakeable, 0 for veto
    int lep4_id;
    int lep4_idx;
    int lep4_mcid;
    int lep4_mcidx;
    int lep4_quality;
    LorentzVector lep4_p4;

    //Isolation
    float lep1_iso;
    float lep2_iso;
    // float lep1_tkIso;
    // float lep2_tkIso;

    int higgs_mass;
    // //Gen Leptons
    int ncharginos;
    int nhiggs;

    vector <LorentzVector> genps_p4;
    vector <LorentzVector> genjets_p4;
    vector <int> genps_id;
    vector <int> genps_id_mother;
    vector <int> genps_idx_mother;
    vector <int> genps_status;
    vector <int> genps_id_grandma;
    vector <LorentzVector> all_leps_p4;
    vector <int> all_leps_id;
    vector <bool> all_leps_veto;
    vector <int> all_leps_mc3idx;

    //Leptons pass numerator ID
    bool lep1_passes_id;
    bool lep2_passes_id;
    bool lep3_passes_id;
    bool lep4_passes_id;

    //Lep1 IDISO
    bool lep1_tight;
    bool lep1_veto;
    bool lep1_fo;

    //Lep2 IDISO
    bool lep2_tight;
    bool lep2_veto;
    bool lep2_fo;

    //Lep3 IDISO
    bool lep3_tight;
    bool lep3_veto;
    bool lep3_fo;

    //Lep4 IDISO
    bool lep4_tight;
    bool lep4_veto;
    bool lep4_fo;

    // bool lep1_passes_MVA;
    // bool lep2_passes_MVA;

    //Imparct parameter
    float lep1_dxyPV;
    float lep2_dxyPV;
    float lep1_d0_err;
    float lep2_d0_err;
    float lep1_dZ;
    float lep2_dZ;
    float lep1_ip3d;
    float lep1_MVA;
    float lep2_MVA;
    // float lep1_MVA_miniaod;
    // float lep2_MVA_miniaod;
    float lep1_ip3d_err;
    float lep2_ip3d;
    float lep2_ip3d_err;

    //nVetoElectrons passing certain pT cuts
    int nVetoElectrons7;
    int nVetoElectrons10;
    int nVetoElectrons25;
    int nVetoMuons5;
    int nVetoMuons10;
    int nVetoMuons25;

    // from miniaod directly
    float lep1_ptrel_ma;
    float lep1_ptratio_ma;
    float lep2_ptrel_ma;
    float lep2_ptratio_ma;

    //PtRel, v0 and v1 (v1 is lepton-subtracted)
    float lep1_ptrel_v1;
    float lep2_ptrel_v1;
    float lep1_ptratio;
    float lep2_ptratio;

    //MiniIso
    float lep1_miniIso;
    float lep2_miniIso;

    //Electron ID variables
    vector <int>   eleID_kfhits;
    vector <float> eleID_oldsigmaietaieta;
    vector <float> eleID_oldsigmaiphiiphi;
    vector <float> eleID_oldcircularity;
    vector <float> eleID_oldr9;
    vector <float> eleID_scletawidth;
    vector <float> eleID_sclphiwidth;
    vector <float> eleID_he;
    vector <float> eleID_psEoverEraw;
    vector <float> eleID_kfchi2;
    vector <float> eleID_chi2_hits;
    vector <float> eleID_fbrem;
    vector <float> eleID_ep;
    vector <float> eleID_eelepout;
    vector <float> eleID_IoEmIop;
    vector <float> eleID_deltaetain;
    vector <float> eleID_deltaphiin;
    vector <float> eleID_deltaetaseed;
    vector <float> eleID_pT;
    vector <bool>  eleID_isbarrel;
    vector <bool>  eleID_isendcap;
    vector <float> eleID_scl_eta;

    //Muon ID variables
    vector <float> muID_dzPV;
    vector <float> muID_ptSig;
    vector <float> muID_ip3dSig;
    vector <float> muID_medMuonPOG;
    vector <float> muID_pt;
    vector <float> muID_eta;

    //InSituFR
    bool lep1_isGoodLeg;
    bool lep2_isGoodLeg;
    bool lep1_isFakeLeg;
    bool lep2_isFakeLeg;
    bool truth_inSituFR;
    bool lep1_multiIso;
    bool lep2_multiIso;
    float lep1_sip;
    float lep2_sip;
    LorentzVector lep1_closeJet;
    LorentzVector lep2_closeJet;
    bool passed_id_inSituFR_lep1;
    bool passed_id_inSituFR_lep2;

    //Trigger matching
    bool fired_trigger;
    bool fired_trigger_ss;
    unsigned int triggers;
    bool lep1_trigMatch_noIsoReq;
    bool lep1_trigMatch_isoReq;
    bool lep2_trigMatch_noIsoReq;
    bool lep2_trigMatch_isoReq;

    //FKW MET
    float met3p0;
    float metphi3p0;

    //MET Filters
    bool passes_met_filters;
    bool passes_chargedcand_filter;
    bool failsRA2Filter;

    // For reminiaod
    float evt_egclean_pfmet;
    float evt_muegclean_pfmet;
    float evt_muegcleanfix_pfmet;
    float evt_uncorr_pfmet;
    bool filt_noBadMuons;
    bool filt_duplicateMuons;
    bool filt_badMuons;

    //Passes 3rd Z/G vetos
    bool madeExtraZ;
    bool madeExtraG;

    //Btag weights
    float weight_btagsf;
    float weight_btagsf_UP;
    float weight_btagsf_DN;
    float weight_btagsf1;
    float weight_btagsf1_UP;
    float weight_btagsf1_DN;
    float weight_btagsf2;
    float weight_btagsf2_UP;
    float weight_btagsf2_DN;
    float weight_btagsf3;
    float weight_btagsf3_UP;
    float weight_btagsf3_DN;
    float weight_btagsf4;
    float weight_btagsf4_UP;
    float weight_btagsf4_DN;

    float weight_btagsf_light_UP;
    float weight_btagsf_light_DN;
    float weight_btagsf1_light_UP;
    float weight_btagsf1_light_DN;
    float weight_btagsf2_light_UP;
    float weight_btagsf2_light_DN;
    float weight_btagsf3_light_UP;
    float weight_btagsf3_light_DN;
    float weight_btagsf4_light_UP;
    float weight_btagsf4_light_DN;
    float weight_btagsf_heavy_UP;
    float weight_btagsf_heavy_DN;
    float weight_btagsf1_heavy_UP;
    float weight_btagsf1_heavy_DN;
    float weight_btagsf2_heavy_UP;
    float weight_btagsf2_heavy_DN;
    float weight_btagsf3_heavy_UP;
    float weight_btagsf3_heavy_DN;
    float weight_btagsf4_heavy_UP;
    float weight_btagsf4_heavy_DN;

    float weight_btagsf_iter_central;
    float weight_btagsf_iter_down_cferr1;
    float weight_btagsf_iter_down_cferr2;
    float weight_btagsf_iter_down_hf;
    float weight_btagsf_iter_down_hfstats1;
    float weight_btagsf_iter_down_hfstats2;
    float weight_btagsf_iter_down_jes;
    float weight_btagsf_iter_down_lf;
    float weight_btagsf_iter_down_lfstats1;
    float weight_btagsf_iter_down_lfstats2;
    float weight_btagsf_iter_up_cferr1;
    float weight_btagsf_iter_up_cferr2;
    float weight_btagsf_iter_up_hf;
    float weight_btagsf_iter_up_hfstats1;
    float weight_btagsf_iter_up_hfstats2;
    float weight_btagsf_iter_up_jes;
    float weight_btagsf_iter_up_lf;
    float weight_btagsf_iter_up_lfstats1;
    float weight_btagsf_iter_up_lfstats2;

    float weight_scale_UP;
    float weight_scale_DN;
    float weight_pdf_UP;
    float weight_pdf_DN;
    float weight_alphas_UP;
    float weight_alphas_DN;
    float weight_isrvar_UP;
    float weight_isrvar_DN;
    float weight_fsrvar_UP;
    float weight_fsrvar_DN;

    // float newweight_scale_UP;
    // float newweight_scale_DN;
    // float newweight_pdf_UP;
    // float newweight_pdf_DN;
    // float newweight_alphas_UP;
    // float newweight_alphas_DN;

    float weight_lepsf1;
    float weight_lepsf2;
    float weight_lepsf3;
    float weight_lepsf;
    float weight_triggersf;
    float weight_pu;
    float weight_isrsf;
    float weight;

    int year;
    float yearlumi;

    //ISR unc #s - 2015
    LorentzVector gl1_p4;
    LorentzVector gl2_p4;
    float glglpt;
    float isr_unc;

    //ISR unc #s - 2016
    int nisrMatch;
    float weight_isr;
    float weight_isr_dy;
    float weight_isr_tt;
    float weight_isr_UP;
    float weight_isr_DN;

    //Lep 3 stuff
    bool lep3_el_etaSC;
    bool lep3_el_conv_vtx_flag;
    bool lep3_el_exp_innerlayers;
    bool lep3_el_threeChargeAgree;
    bool lep3_el_dxyPV;
    bool lep3_el_ip3d;
    bool lep3_el_dzPV;
    float lep3_el_MVA_value;
    bool lep3_el_MVA;
    bool lep3_iso_RA5;
    bool lep3_iso_RA7;
    bool lep3_passes_RA5;
    bool lep3_passes_RA7;
    bool lep3_passes_RA7v2;
    bool lep3_mu_dxyPV;
    bool lep3_mu_ip3d;
    bool lep3_mu_dzPV;
    bool lep3_mu_ptErr;
    bool lep3_isTrigSafeNoIsov1;
    bool lep3_isTrigSafev1;

    //Lep 4 stuff
    bool lep4_el_etaSC;
    bool lep4_el_conv_vtx_flag;
    bool lep4_el_exp_innerlayers;
    bool lep4_el_threeChargeAgree;
    bool lep4_el_dxyPV;
    bool lep4_el_ip3d;
    bool lep4_el_dzPV;
    float lep4_el_MVA_value;
    bool lep4_el_MVA;
    bool lep4_iso_RA5;
    bool lep4_iso_RA7;
    bool lep4_passes_RA5;
    bool lep4_passes_RA7;
    bool lep4_mu_dxyPV;
    bool lep4_mu_ip3d;
    bool lep4_mu_dzPV;
    bool lep4_mu_ptErr;
    bool lep4_isTrigSafeNoIsov1;
    bool lep4_isTrigSafev1;

    // trijet stuff
    int ntrijets;
    vector<float> trijet_discs;
    int trijet_njetsnonb;
    int trijet_njetsb;
    float trijet_meandisc;
    float trijet_leadingdisc;
    float trijet_subleadingdisc;
    int trijet_numhigh;
    float trijet_frachigh;

    // maybe bdt variables?
    int bdt_nforwardjets20;
    float bdt_avgcdisc;


    // bdt variables
    float bdt_nbtags;
    float bdt_njets;
    float bdt_met;
    float bdt_ptl2;
    float bdt_nlb40;
    float bdt_ntb40;
    float bdt_nleps;
    float bdt_htb;
    float bdt_ml1j1;
    float bdt_dphil1l2;
    float bdt_maxmjoverpt;
    float bdt_detal1l2;
    float bdt_q1;
    float bdt_ptj1;
    float bdt_ptj6;
    float bdt_ptj7;
    float bdt_ptj8;
    float bdt_ptl1;
    float bdt_ptl3;

    float bdt_jec_up_nbtags;
    float bdt_jec_dn_nbtags;
    float bdt_jer_up_nbtags;
    float bdt_jer_dn_nbtags;

    float bdt_jec_up_njets;
    float bdt_jec_dn_njets;
    float bdt_jer_up_njets;
    float bdt_jer_dn_njets;

    float bdt_jec_up_met;
    float bdt_jec_dn_met;
    float bdt_jer_up_met;
    float bdt_jer_dn_met;

    float bdt_jec_up_htb;
    float bdt_jec_dn_htb;
    float bdt_jer_up_htb;
    float bdt_jer_dn_htb;

    float bdt_jec_up_nlb40;
    float bdt_jec_dn_nlb40;
    float bdt_jer_up_nlb40;
    float bdt_jer_dn_nlb40;

    float bdt_jec_up_ntb40;
    float bdt_jec_dn_ntb40;
    float bdt_jer_up_ntb40;
    float bdt_jer_dn_ntb40;

    float bdt_disc;
    float bdt_disc_jec_up;
    float bdt_disc_jer_up;
    float bdt_disc_jec_dn;
    float bdt_disc_jer_dn;

    int bdt_run2_SRDISC;
    float bdt_run2_disc;
    float bdt_run2_disc_jec_up;
    float bdt_run2_disc_jer_up;
    float bdt_run2_disc_jec_dn;
    float bdt_run2_disc_jer_dn;

};

#endif
