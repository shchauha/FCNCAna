#include "helper_babymaker.h"
#include "CORE/Tools/btagsf/BTagCalibrationStandalone.h"
#include "CORE/Tools/jetcorr/Utilities.icc"
#include "CORE/Tools/jetcorr/SimpleJetCorrectionUncertainty.h"
#include "CORE/Tools/jetcorr/JetCorrectionUncertainty.h"

#include "misc/signal_regions.h"
#include "misc/bdt.h"
namespace bdtrun2 {
#include "misc/bdt_run2.h"
}
#include "misc/common_utils.h"
// #include "misc/toptagger.h"

using namespace tas;

const bool applyBtagSFs = true;

//Main functions
void babyMaker::MakeBabyNtuple(const char* output_name, bool isFastsim, int iSignal){

  //Create Baby
  BabyFile = new TFile(Form("%s/%s", path.Data(), output_name), "RECREATE");
  BabyFile->cd();
  BabyTree = new TTree("t", "SS2015 Baby Ntuple");

  //Define Branches
  BabyTree->Branch("lep1_isPrompt"                                           , &lep1_isPrompt                                           );
  BabyTree->Branch("lep1_isDirectPrompt"                                     , &lep1_isDirectPrompt                                     );
  BabyTree->Branch("lep1_isStat3"                                            , &lep1_isStat3                                            );
  BabyTree->Branch("lep2_isPrompt"                                           , &lep2_isPrompt                                           );
  BabyTree->Branch("lep2_isDirectPrompt"                                     , &lep2_isDirectPrompt                                     );
  BabyTree->Branch("lep2_isStat3"                                            , &lep2_isStat3                                            );
  BabyTree->Branch("met"                                                     , &met                                                     );
  BabyTree->Branch("metPhi"                                                  , &metPhi                                                  );
  BabyTree->Branch("rawmet"                                                  , &rawmet                                                  );
  BabyTree->Branch("rawmetPhi"                                               , &rawmetPhi                                               );
  BabyTree->Branch("calomet"                                                  , &calomet                                                  );
  BabyTree->Branch("calometPhi"                                               , &calometPhi                                               );
  BabyTree->Branch("event"                                                   , &event                                                   );
  BabyTree->Branch("lumi"                                                    , &lumi                                                    );
  BabyTree->Branch("run"                                                     , &run                                                     );
  BabyTree->Branch("filt_csc"                                                , &filt_csc                                                );
  BabyTree->Branch("filt_hbhe"                                               , &filt_hbhe                                               );
  BabyTree->Branch("filt_hcallaser"                                          , &filt_hcallaser                                          );
  BabyTree->Branch("filt_ecaltp"                                             , &filt_ecaltp                                             );
  BabyTree->Branch("filt_trkfail"                                            , &filt_trkfail                                            );
  BabyTree->Branch("is_real_data"                                            , &is_real_data                                            );
  BabyTree->Branch("scale1fb"                                                , &scale1fb                                                );
  BabyTree->Branch("genps_weight"                                                , &genps_weight                                                );
  BabyTree->Branch("xsec"                                                    , &xsec                                                    );
  BabyTree->Branch("neventstotal"                                                    , &neventstotal                                                    );
  BabyTree->Branch("xsec_ps"                                                    , &xsec_ps                                                    );
  BabyTree->Branch("sparmNames"                                              , &sparmNames                                              );
  BabyTree->Branch("sparms"                                                  , &sparms                                                  );
  BabyTree->Branch("xsec_error"                                              , &xsec_error                                              );
  BabyTree->Branch("kfactor"                                                 , &kfactor                                                 );
  BabyTree->Branch("gen_met"                                                 , &gen_met                                                 );
  // BabyTree->Branch("genweights"                                              , &genweights                                              );
  BabyTree->Branch("pdfweights"                                              , &pdfweights                                              );
  // BabyTree->Branch("genweightsID"                                            , &genweightsID                                            );
  BabyTree->Branch("gen_met_phi"                                             , &gen_met_phi                                             );
  BabyTree->Branch("skim"                                                   , &skim                                                   );
  BabyTree->Branch("skim_nomet"                                                   , &skim_nomet                                                   );
  BabyTree->Branch("nleps"                                                   , &nleps                                                   );
  BabyTree->Branch("sr"                                                   , &sr                                                   );
  BabyTree->Branch("br"                                                   , &br                                                   );
  BabyTree->Branch("ss_nleps"                                                   , &ss_nleps                                                   );
  BabyTree->Branch("ss_br"                                                   , &ss_br                                                   );
  BabyTree->Branch("ss_sr"                                                   , &ss_sr                                                   );
  BabyTree->Branch("ss_region"                                                   , &ss_region                                                   );
  BabyTree->Branch("njets"                                                   , &njets                                                   );
  BabyTree->Branch("njetsAG"                                                   , &njetsAG                                                   );
  BabyTree->Branch("nbtagsAG"                                                   , &nbtagsAG                                                   );
  BabyTree->Branch("njets_raw"                                               , &njets_raw                                               );
  BabyTree->Branch("hyp_class"                                               , &hyp_class                                               );
  BabyTree->Branch("lep1_p4"                                                 , &lep1_p4                                                 );
  BabyTree->Branch("lep2_p4"                                                 , &lep2_p4                                                 );
  BabyTree->Branch("lep1_pt"                                                 , &lep1_pt                                                 );
  BabyTree->Branch("lep2_pt"                                                 , &lep2_pt                                                 );
  BabyTree->Branch("lep3_pt"                                                 , &lep3_pt                                                 );
  BabyTree->Branch("lep1_eta"                                                 , &lep1_eta                                                 );
  BabyTree->Branch("lep2_eta"                                                 , &lep2_eta                                                 );
  BabyTree->Branch("lep3_eta"                                                 , &lep3_eta                                                 );
  BabyTree->Branch("lep1_phi"                                                 , &lep1_phi                                                 );
  BabyTree->Branch("lep2_phi"                                                 , &lep2_phi                                                 );
  BabyTree->Branch("lep3_phi"                                                 , &lep3_phi                                                 );
  BabyTree->Branch("ht_raw"                                                  , &ht_raw                                                  );
  BabyTree->Branch("ht"                                                      , &ht                                                      );
  BabyTree->Branch("lep1_motherID"                                           , &lep1_motherID                                           );
  BabyTree->Branch("lep2_motherID"                                           , &lep2_motherID                                           );
  BabyTree->Branch("lep1_mc_id"                                              , &lep1_mc_id                                              );
  BabyTree->Branch("lep2_mc_id"                                              , &lep2_mc_id                                              );
  BabyTree->Branch("lep1_mc_motherid"                                              , &lep1_mc_motherid                                              );
  BabyTree->Branch("lep2_mc_motherid"                                              , &lep2_mc_motherid                                              );
  BabyTree->Branch("lep1_id"                                                 , &lep1_id                                                 );
  BabyTree->Branch("lep2_id"                                                 , &lep2_id                                                 );
  BabyTree->Branch("lep1_coneCorrPt"                                         , &lep1_coneCorrPt                                         );
  BabyTree->Branch("lep2_coneCorrPt"                                         , &lep2_coneCorrPt                                         );
  BabyTree->Branch("lep3_coneCorrPt"                                         , &lep3_coneCorrPt                                         );
  BabyTree->Branch("lep4_coneCorrPt"                                         , &lep4_coneCorrPt                                         );
  BabyTree->Branch("lep1_idx"                                                , &lep1_idx                                                );
  BabyTree->Branch("lep2_idx"                                                , &lep2_idx                                                );
  BabyTree->Branch("jets"                                                    , &jets                                                    );
  BabyTree->Branch("jets_jec_up"                                            , &jets_jec_up                                            );
  BabyTree->Branch("jets_jec_dn"                                            , &jets_jec_dn                                            );
  BabyTree->Branch("jets_disc_up"                                            , &jets_disc_up                                            );
  BabyTree->Branch("jets_disc_dn"                                            , &jets_disc_dn                                            );
  BabyTree->Branch("bjets_jec_up"                                            , &bjets_jec_up                                            );
  BabyTree->Branch("bjets_jec_dn"                                            , &bjets_jec_dn                                            );
  BabyTree->Branch("btags_flavor"                                              , &btags_flavor                                              );
  BabyTree->Branch("btags_disc"                                              , &btags_disc                                              );
  BabyTree->Branch("btags_cdisc"                                              , &btags_cdisc                                              );
  BabyTree->Branch("jets_flavor"                                               , &jets_flavor                                               );
  BabyTree->Branch("jets_disc"                                               , &jets_disc                                               );
  BabyTree->Branch("jets_bsf"                                               , &jets_bsf                                               );
  BabyTree->Branch("jets_JEC"                                                , &jets_JEC                                                );
  BabyTree->Branch("btags_JEC"                                               , &btags_JEC                                               );
  BabyTree->Branch("jets_undoJEC"                                            , &jets_undoJEC                                            );
  BabyTree->Branch("btags_undoJEC"                                           , &btags_undoJEC                                           );
  BabyTree->Branch("btags_unc"                                               , &btags_unc                                               );
  BabyTree->Branch("btags_eff"                                               , &btags_eff                                               );
  BabyTree->Branch("btags_effpt"                                               , &btags_effpt                                               );
  BabyTree->Branch("btags_sf"                                               , &btags_sf                                               );
  BabyTree->Branch("jets_unc"                                                , &jets_unc                                                );
  BabyTree->Branch("btags"                                                   , &btags                                                   );
  BabyTree->Branch("nbtags"                                                  , &nbtags                                                  );
  BabyTree->Branch("nbtags_raw"                                              , &nbtags_raw                                              );
  BabyTree->Branch("sf_dilepTrig_hpt"                                        , &sf_dilepTrig_hpt                                        );
  BabyTree->Branch("sf_dilepTrig_lpt"                                        , &sf_dilepTrig_lpt                                        );
  BabyTree->Branch("sf_dilepTrig_vlpt"                                       , &sf_dilepTrig_vlpt                                       );
  BabyTree->Branch("hyp_type"                                                , &hyp_type                                                );
  BabyTree->Branch("sf_dilep_eff"                                            , &sf_dilep_eff                                            );
  BabyTree->Branch("mt"                                                      , &mt                                                      );
  BabyTree->Branch("mt_l2"                                                   , &mt_l2                                                   );
  BabyTree->Branch("mt2"                                                     , &mt2                                                     );
  BabyTree->Branch("mtmin"                                                   , &mtmin                                                   );
  BabyTree->Branch("lep1_id_gen"                                             , &lep1_id_gen                                             );
  BabyTree->Branch("lep2_id_gen"                                             , &lep2_id_gen                                             );
  BabyTree->Branch("lep1_p4_gen"                                             , &lep1_p4_gen                                             );
  BabyTree->Branch("lep2_p4_gen"                                             , &lep2_p4_gen                                             );
  BabyTree->Branch("lep3_id"                                                 , &lep3_id                                                 );
  BabyTree->Branch("lep3_idx"                                                , &lep3_idx                                                );
  BabyTree->Branch("lep3_p4"                                                 , &lep3_p4                                                 );
  BabyTree->Branch("lep3_quality"                                            , &lep3_quality                                            );
  BabyTree->Branch("lep4_id"                                                 , &lep4_id                                                 );
  BabyTree->Branch("lep4_idx"                                                , &lep4_idx                                                );
  BabyTree->Branch("lep4_p4"                                                 , &lep4_p4                                                 );
  BabyTree->Branch("lep1_iso"                                                , &lep1_iso                                                );
  BabyTree->Branch("lep2_iso"                                                , &lep2_iso                                                );
  // BabyTree->Branch("lep1_tkIso"                                              , &lep1_tkIso                                              );
  // BabyTree->Branch("lep2_tkIso"                                              , &lep2_tkIso                                              );
  BabyTree->Branch("dilep_p4"                                                , &dilep_p4                                                );
  BabyTree->Branch("dilep_mass"                                                , &dilep_mass                                                );
  BabyTree->Branch("mass13"                                                , &mass13                                                );
  BabyTree->Branch("mass23"                                                , &mass23                                                );
  BabyTree->Branch("mass123"                                                , &mass123                                                );
  BabyTree->Branch("ncharginos"                                                , &ncharginos                                                );
  BabyTree->Branch("nhiggs"                                                , &nhiggs                                                );
  BabyTree->Branch("higgs_mass"                                                , &higgs_mass                                                );
  BabyTree->Branch("genps_p4"                                                , &genps_p4                                                );
  BabyTree->Branch("genjets_p4"                                                , &genjets_p4                                                );
  BabyTree->Branch("genps_id"                                                , &genps_id                                                );
  BabyTree->Branch("genps_id_mother"                                         , &genps_id_mother                                         );
  BabyTree->Branch("genps_idx_mother"                                        , &genps_idx_mother                                        );
  BabyTree->Branch("genps_status"                                            , &genps_status                                            );
  BabyTree->Branch("genps_id_grandma"                                        , &genps_id_grandma                                        );

  BabyTree->Branch("all_leps_p4"                                        , &all_leps_p4                                        );
  BabyTree->Branch("all_leps_id"                                        , &all_leps_id                                        );
  BabyTree->Branch("all_leps_veto"                                        , &all_leps_veto                                        );
  BabyTree->Branch("all_leps_mc3idx"                                        , &all_leps_mc3idx                                        );

  BabyTree->Branch("lep1_passes_id"                                          , &lep1_passes_id                                          );
  BabyTree->Branch("lep2_passes_id"                                          , &lep2_passes_id                                          );
  BabyTree->Branch("lep3_passes_id"                                          , &lep3_passes_id                                          );
  BabyTree->Branch("lep4_passes_id"                                          , &lep4_passes_id                                          );
  BabyTree->Branch("lep1_tight"                                              , &lep1_tight                                              );
  BabyTree->Branch("lep1_veto"                                               , &lep1_veto                                               );
  BabyTree->Branch("lep1_fo"                                                 , &lep1_fo                                                 );
  BabyTree->Branch("lep2_tight"                                              , &lep2_tight                                              );
  BabyTree->Branch("lep2_veto"                                               , &lep2_veto                                               );
  BabyTree->Branch("lep2_fo"                                                 , &lep2_fo                                                 );
  BabyTree->Branch("lep3_tight"                                              , &lep3_tight                                              );
  BabyTree->Branch("lep3_veto"                                               , &lep3_veto                                               );
  BabyTree->Branch("lep3_fo"                                                 , &lep3_fo                                                 );
  BabyTree->Branch("lep4_tight"                                              , &lep4_tight                                              );
  BabyTree->Branch("lep4_veto"                                               , &lep4_veto                                               );
  BabyTree->Branch("lep4_fo"                                                 , &lep4_fo                                                 );
  BabyTree->Branch("lep1_dxyPV"                                              , &lep1_dxyPV                                              );
  BabyTree->Branch("lep2_dxyPV"                                              , &lep2_dxyPV                                              );
  BabyTree->Branch("lep1_dZ"                                                 , &lep1_dZ                                                 );
  BabyTree->Branch("lep2_dZ"                                                 , &lep2_dZ                                                 );
  BabyTree->Branch("lep1_d0_err"                                             , &lep1_d0_err                                             );
  BabyTree->Branch("lep2_d0_err"                                             , &lep2_d0_err                                             );
  BabyTree->Branch("lep1_ip3d"                                               , &lep1_ip3d                                               );
  BabyTree->Branch("lep2_ip3d"                                               , &lep2_ip3d                                               );
  BabyTree->Branch("lep1_MVA"                                                , &lep1_MVA                                                );
  BabyTree->Branch("lep2_MVA"                                                , &lep2_MVA                                                );
  // BabyTree->Branch("lep1_MVA_miniaod"                                        , &lep1_MVA_miniaod                                        );
  // BabyTree->Branch("lep2_MVA_miniaod"                                        , &lep2_MVA_miniaod                                        );
  // BabyTree->Branch("lep1_passes_MVA"                                        , &lep1_passes_MVA                                        );
  // BabyTree->Branch("lep2_passes_MVA"                                        , &lep2_passes_MVA                                        );
  BabyTree->Branch("lep1_ip3d_err"                                           , &lep1_ip3d_err                                           );
  BabyTree->Branch("lep2_ip3d_err"                                           , &lep2_ip3d_err                                           );
  BabyTree->Branch("nVetoElectrons7"                                         , &nVetoElectrons7                                         );
  BabyTree->Branch("nVetoElectrons10"                                        , &nVetoElectrons10                                        );
  BabyTree->Branch("nVetoElectrons25"                                        , &nVetoElectrons25                                        );
  BabyTree->Branch("nVetoMuons5"                                             , &nVetoMuons5                                             );
  BabyTree->Branch("nVetoMuons10"                                            , &nVetoMuons10                                            );
  BabyTree->Branch("nVetoMuons25"                                            , &nVetoMuons25                                            );
  BabyTree->Branch("filename"                                                , &filename                                                );
  BabyTree->Branch("lep1_ptrel_ma"                                           , &lep1_ptrel_ma                                           );
  BabyTree->Branch("lep1_ptratio_ma"                                           , &lep1_ptratio_ma                                           );
  BabyTree->Branch("lep2_ptrel_ma"                                           , &lep2_ptrel_ma                                           );
  BabyTree->Branch("lep2_ptratio_ma"                                           , &lep2_ptratio_ma                                           );
  BabyTree->Branch("lep1_ptrel_v1"                                           , &lep1_ptrel_v1                                           );
  BabyTree->Branch("lep2_ptrel_v1"                                           , &lep2_ptrel_v1                                           );
  BabyTree->Branch("lep1_ptratio"                                           , &lep1_ptratio                                           );
  BabyTree->Branch("lep2_ptratio"                                           , &lep2_ptratio                                           );
  BabyTree->Branch("lep1_miniIso"                                            , &lep1_miniIso                                            );
  BabyTree->Branch("lep2_miniIso"                                            , &lep2_miniIso                                            );
  BabyTree->Branch("trueNumInt"                                              , &trueNumInt                                              );
  BabyTree->Branch("nGoodVertices"                                           , &nGoodVertices                                           );
  BabyTree->Branch("lep1_trigMatch_noIsoReq"                                 , &lep1_trigMatch_noIsoReq                                 );
  BabyTree->Branch("lep1_trigMatch_isoReq"                                   , &lep1_trigMatch_isoReq                                   );
  BabyTree->Branch("lep2_trigMatch_noIsoReq"                                 , &lep2_trigMatch_noIsoReq                                 );
  BabyTree->Branch("lep2_trigMatch_isoReq"                                   , &lep2_trigMatch_isoReq                                   );
  BabyTree->Branch("passes_met_filters"                                      , &passes_met_filters                                      );
  BabyTree->Branch("passes_chargedcand_filter"                                      , &passes_chargedcand_filter                                      );
  BabyTree->Branch("evt_egclean_pfmet"                                      , &evt_egclean_pfmet                                      );
  BabyTree->Branch("evt_muegclean_pfmet"                                      , &evt_muegclean_pfmet                                      );
  BabyTree->Branch("evt_muegcleanfix_pfmet"                                      , &evt_muegcleanfix_pfmet                                      );
  BabyTree->Branch("evt_uncorr_pfmet"                                      , &evt_uncorr_pfmet                                      );
  BabyTree->Branch("filt_noBadMuons"                                      , &filt_noBadMuons                                      );
  BabyTree->Branch("filt_duplicateMuons"                                      , &filt_duplicateMuons                                      );
  BabyTree->Branch("filt_badMuons"                                      , &filt_badMuons                                      );
  BabyTree->Branch("failsRA2Filter"                                      , &failsRA2Filter                                      );
  BabyTree->Branch("madeExtraZ"                                              , &madeExtraZ                                              );
  BabyTree->Branch("madeExtraG"                                              , &madeExtraG                                              );
  BabyTree->Branch("lep3_mcid"                                               , &lep3_mcid                                               );
  BabyTree->Branch("lep3_mc_motherid"                                               , &lep3_mc_motherid                                               );
  BabyTree->Branch("lep3_mc3idx"                                              , &lep3_mc3idx                                              );
  BabyTree->Branch("lep3_motherID"                                              , &lep3_motherID                                              );
  BabyTree->Branch("lep4_mcid"                                               , &lep4_mcid                                               );
  BabyTree->Branch("lep4_mcidx"                                              , &lep4_mcidx                                              );

  BabyTree->Branch("njets_unc_up"                                            , &njets_unc_up                                            );
  BabyTree->Branch("njets_unc_dn"                                            , &njets_unc_dn                                            );
  BabyTree->Branch("ht_unc_up"                                               , &ht_unc_up                                               );
  BabyTree->Branch("ht_unc_dn"                                               , &ht_unc_dn                                               );
  BabyTree->Branch("nbtags_unc_up"                                           , &nbtags_unc_up                                           );
  BabyTree->Branch("nbtags_unc_dn"                                           , &nbtags_unc_dn                                           );
  BabyTree->Branch("met_unc_up"                                              , &met_unc_up                                              );
  BabyTree->Branch("met_unc_dn"                                              , &met_unc_dn                                              );
  BabyTree->Branch("metPhi_unc_up"                                           , &metPhi_unc_up                                           );
  BabyTree->Branch("metPhi_unc_dn"                                           , &metPhi_unc_dn                                           );
  BabyTree->Branch("njets_JER_up"                                            , &njets_JER_up                                            );
  BabyTree->Branch("njets_JER_dn"                                            , &njets_JER_dn                                            );
  BabyTree->Branch("ht_JER_up"                                               , &ht_JER_up                                               );
  BabyTree->Branch("ht_JER_dn"                                               , &ht_JER_dn                                               );
  BabyTree->Branch("nbtags_JER_up"                                           , &nbtags_JER_up                                           );
  BabyTree->Branch("nbtags_JER_dn"                                           , &nbtags_JER_dn                                           );
  BabyTree->Branch("met_JER_up"                                              , &met_JER_up                                              );
  BabyTree->Branch("met_JER_dn"                                              , &met_JER_dn                                              );
  BabyTree->Branch("metPhi_JER_up"                                           , &metPhi_JER_up                                           );
  BabyTree->Branch("metPhi_JER_dn"                                           , &metPhi_JER_dn                                           );
  BabyTree->Branch("lep2_genps_isHardProcess"                                , &lep2_genps_isHardProcess                                );
  BabyTree->Branch("lep2_genps_fromHardProcessFinalState"                    , &lep2_genps_fromHardProcessFinalState                    );
  BabyTree->Branch("lep2_genps_fromHardProcessDecayed"                       , &lep2_genps_fromHardProcessDecayed                       );
  BabyTree->Branch("lep2_genps_isDirectHardProcessTauDecayProductFinalState" , &lep2_genps_isDirectHardProcessTauDecayProductFinalState );
  BabyTree->Branch("lep2_genps_fromHardProcessBeforeFSR"                     , &lep2_genps_fromHardProcessBeforeFSR                     );
  BabyTree->Branch("lep2_genps_isLastCopy"                                   , &lep2_genps_isLastCopy                                   );
  BabyTree->Branch("lep2_genps_isLastCopyBeforeFSR"                          , &lep2_genps_isLastCopyBeforeFSR                          );
  BabyTree->Branch("lep1_genps_isHardProcess"                                , &lep1_genps_isHardProcess                                );
  BabyTree->Branch("lep1_genps_fromHardProcessFinalState"                    , &lep1_genps_fromHardProcessFinalState                    );
  BabyTree->Branch("lep1_genps_fromHardProcessDecayed"                       , &lep1_genps_fromHardProcessDecayed                       );
  BabyTree->Branch("lep1_genps_isDirectHardProcessTauDecayProductFinalState" , &lep1_genps_isDirectHardProcessTauDecayProductFinalState );
  BabyTree->Branch("lep1_genps_fromHardProcessBeforeFSR"                     , &lep1_genps_fromHardProcessBeforeFSR                     );
  BabyTree->Branch("lep1_genps_isLastCopy"                                   , &lep1_genps_isLastCopy                                   );
  BabyTree->Branch("lep1_genps_isLastCopyBeforeFSR"                          , &lep1_genps_isLastCopyBeforeFSR                          );
  BabyTree->Branch("lep1_mc3idx"                                             , &lep1_mc3idx                                             );
  BabyTree->Branch("lep2_mc3idx"                                             , &lep2_mc3idx                                             );
  BabyTree->Branch("lep1_el_conv_vtx_flag"   , &lep1_el_conv_vtx_flag   );
  BabyTree->Branch("lep2_el_conv_vtx_flag"   , &lep2_el_conv_vtx_flag   );
  BabyTree->Branch("lep1_el_exp_innerlayers" , &lep1_el_exp_innerlayers );
  BabyTree->Branch("lep1_mu_ptErr" , &lep1_mu_ptErr );
  BabyTree->Branch("lep2_mu_ptErr" , &lep2_mu_ptErr );
  BabyTree->Branch("lep1_el_threeChargeAgree" , &lep1_el_threeChargeAgree );
  BabyTree->Branch("lep2_el_threeChargeAgree" , &lep2_el_threeChargeAgree );
  BabyTree->Branch("lep2_el_exp_innerlayers" , &lep2_el_exp_innerlayers );
  BabyTree->Branch("lep1_el_etaSC" , &lep1_el_etaSC );
  BabyTree->Branch("lep2_el_etaSC" , &lep2_el_etaSC );
  BabyTree->Branch("lep1_nPixelMiss"   , &lep1_nPixelMiss   );
  BabyTree->Branch("lep2_nPixelMiss"   , &lep2_nPixelMiss   );
  BabyTree->Branch("lep1_tightCharge"   , &lep1_tightCharge   );
  BabyTree->Branch("lep2_tightCharge"   , &lep2_tightCharge   );
  BabyTree->Branch("lep1_isTrigSafeNoIsov1"  , &lep1_isTrigSafeNoIsov1  );
  BabyTree->Branch("lep1_isTrigSafev1"       , &lep1_isTrigSafev1       );
  BabyTree->Branch("lep2_isTrigSafeNoIsov1"  , &lep2_isTrigSafeNoIsov1  );
  BabyTree->Branch("lep2_isTrigSafev1"       , &lep2_isTrigSafev1       );
  BabyTree->Branch("lep1_tightMuonPOG"  , &lep1_tightMuonPOG  );
  BabyTree->Branch("lep1_mediumMuonPOG"  , &lep1_mediumMuonPOG  );
  BabyTree->Branch("lep1_looseMuonPOG"  , &lep1_looseMuonPOG  );
  BabyTree->Branch("lep2_tightMuonPOG"  , &lep2_tightMuonPOG  );
  BabyTree->Branch("lep2_mediumMuonPOG"  , &lep2_mediumMuonPOG  );
  BabyTree->Branch("lep2_looseMuonPOG"  , &lep2_looseMuonPOG  );
  BabyTree->Branch("lep3_mediumMuonPOG"  , &lep3_mediumMuonPOG  );

    // FG: add VVV lepton info
  BabyTree->Branch("lep1_el_vvv_ss_fo_id",&lep1_el_vvv_ss_fo_id);
  BabyTree->Branch("lep1_el_vvv_3l_fo_id",&lep1_el_vvv_3l_fo_id);
  BabyTree->Branch("lep1_el_vvv_ss_id",&lep1_el_vvv_ss_id);
  BabyTree->Branch("lep1_el_vvv_3l_id",&lep1_el_vvv_3l_id);
  BabyTree->Branch("lep1_mu_vvv_fo_id",&lep1_mu_vvv_fo_id);
  BabyTree->Branch("lep1_mu_vvv_id",&lep1_mu_vvv_id);
  BabyTree->Branch("lep1_vvv_iso",&lep1_vvv_iso);
  BabyTree->Branch("lep2_el_vvv_ss_fo_id",&lep2_el_vvv_ss_fo_id);
  BabyTree->Branch("lep2_el_vvv_3l_fo_id",&lep2_el_vvv_3l_fo_id);
  BabyTree->Branch("lep2_el_vvv_ss_id",&lep2_el_vvv_ss_id);
  BabyTree->Branch("lep2_el_vvv_3l_id",&lep2_el_vvv_3l_id);
  BabyTree->Branch("lep2_mu_vvv_fo_id",&lep2_mu_vvv_fo_id);
  BabyTree->Branch("lep2_mu_vvv_id",&lep2_mu_vvv_id);
  BabyTree->Branch("lep2_vvv_iso",&lep2_vvv_iso);

  BabyTree->Branch("prefire_sf"  , &prefire_sf  );
  BabyTree->Branch("prefire_sfup"  , &prefire_sfup  );
  BabyTree->Branch("prefire_sfdown"  , &prefire_sfdown  );
  BabyTree->Branch("prefire2016_sf"  , &prefire2016_sf  );
  BabyTree->Branch("prefire2016_sfup"  , &prefire2016_sfup  );
  BabyTree->Branch("prefire2016_sfdown"  , &prefire2016_sfdown  );
  BabyTree->Branch("prefire2017_sf"  , &prefire2017_sf  );
  BabyTree->Branch("prefire2017_sfup"  , &prefire2017_sfup  );
  BabyTree->Branch("prefire2017_sfdown"  , &prefire2017_sfdown  );
  BabyTree->Branch("prefire2017ele_sf"  , &prefire2017ele_sf  );
  BabyTree->Branch("prefire2017ele_sfup"  , &prefire2017ele_sfup  );
  BabyTree->Branch("prefire2017ele_sfdown"  , &prefire2017ele_sfdown  );


  BabyTree->Branch("is_fastsim", &is_fastsim);

  //InSituFR
  BabyTree->Branch("lep1_isGoodLeg"                                          , &lep1_isGoodLeg                                                                          );
  BabyTree->Branch("lep2_isGoodLeg"                                          , &lep2_isGoodLeg                                                                          );
  BabyTree->Branch("lep1_isFakeLeg"                                          , &lep1_isFakeLeg                                                                          );
  BabyTree->Branch("lep2_isFakeLeg"                                          , &lep2_isFakeLeg                                                                          );
  BabyTree->Branch("lep1_multiIso"                                           , &lep1_multiIso                                                                           );
  BabyTree->Branch("lep2_multiIso"                                           , &lep2_multiIso                                                                           );
  BabyTree->Branch("lep1_sip"                                                , &lep1_sip                                                                                );
  BabyTree->Branch("lep2_sip"                                                , &lep2_sip                                                                                );
  BabyTree->Branch("lep1_closeJet"                                           , &lep1_closeJet                                                                           );
  BabyTree->Branch("lep2_closeJet"                                           , &lep2_closeJet                                                                           );
  BabyTree->Branch("passed_id_inSituFR_lep1"                                 , &passed_id_inSituFR_lep1                                                                 );
  BabyTree->Branch("passed_id_inSituFR_lep2"                                 , &passed_id_inSituFR_lep2                                                                 );

  //Triggers
  BabyTree->Branch("fired_trigger"                                           , &fired_trigger                                                                           );
  BabyTree->Branch("fired_trigger_ss"                                           , &fired_trigger_ss                                                                           );
  BabyTree->Branch("triggers"                                                , &triggers                                                                                );

  BabyTree->Branch("weight_btagsf"                                           , &weight_btagsf                                                                           );
  // BabyTree->Branch("weight_btagsf_UP"                                        , &weight_btagsf_UP                                                                        );
  // BabyTree->Branch("weight_btagsf_DN"                                        , &weight_btagsf_DN                                                                        );
  // BabyTree->Branch("weight_btagsf1"                                           , &weight_btagsf1                                                                           );
  // BabyTree->Branch("weight_btagsf1_UP"                                        , &weight_btagsf1_UP                                                                        );
  // BabyTree->Branch("weight_btagsf1_DN"                                        , &weight_btagsf1_DN                                                                        );
  // BabyTree->Branch("weight_btagsf2"                                           , &weight_btagsf2                                                                           );
  // BabyTree->Branch("weight_btagsf2_UP"                                        , &weight_btagsf2_UP                                                                        );
  // BabyTree->Branch("weight_btagsf2_DN"                                        , &weight_btagsf2_DN                                                                        );
  // BabyTree->Branch("weight_btagsf3"                                           , &weight_btagsf3                                                                           );
  // BabyTree->Branch("weight_btagsf3_UP"                                        , &weight_btagsf3_UP                                                                        );
  // BabyTree->Branch("weight_btagsf3_DN"                                        , &weight_btagsf3_DN                                                                        );
  // BabyTree->Branch("weight_btagsf4"                                           , &weight_btagsf4                                                                           );
  // BabyTree->Branch("weight_btagsf4_UP"                                        , &weight_btagsf4_UP                                                                        );
  // BabyTree->Branch("weight_btagsf4_DN"                                        , &weight_btagsf4_DN                                                                        );
  BabyTree->Branch("weight_btagsf_light_UP"                                        , &weight_btagsf_light_UP                                                                        );
  BabyTree->Branch("weight_btagsf_light_DN"                                        , &weight_btagsf_light_DN                                                                        );
  BabyTree->Branch("weight_btagsf_heavy_UP"                                        , &weight_btagsf_heavy_UP                                                                        );
  BabyTree->Branch("weight_btagsf_heavy_DN"                                        , &weight_btagsf_heavy_DN                                                                        );
  // BabyTree->Branch("weight_btagsf1_light_UP"                                        , &weight_btagsf1_light_UP                                                                        );
  // BabyTree->Branch("weight_btagsf1_light_DN"                                        , &weight_btagsf1_light_DN                                                                        );
  // BabyTree->Branch("weight_btagsf2_light_UP"                                        , &weight_btagsf2_light_UP                                                                        );
  // BabyTree->Branch("weight_btagsf2_light_DN"                                        , &weight_btagsf2_light_DN                                                                        );
  // BabyTree->Branch("weight_btagsf3_light_UP"                                        , &weight_btagsf3_light_UP                                                                        );
  // BabyTree->Branch("weight_btagsf3_light_DN"                                        , &weight_btagsf3_light_DN                                                                        );
  // BabyTree->Branch("weight_btagsf4_light_UP"                                        , &weight_btagsf4_light_UP                                                                        );
  // BabyTree->Branch("weight_btagsf4_light_DN"                                        , &weight_btagsf4_light_DN                                                                        );
  // BabyTree->Branch("weight_btagsf1_heavy_UP"                                        , &weight_btagsf1_heavy_UP                                                                        );
  // BabyTree->Branch("weight_btagsf1_heavy_DN"                                        , &weight_btagsf1_heavy_DN                                                                        );
  // BabyTree->Branch("weight_btagsf2_heavy_UP"                                        , &weight_btagsf2_heavy_UP                                                                        );
  // BabyTree->Branch("weight_btagsf2_heavy_DN"                                        , &weight_btagsf2_heavy_DN                                                                        );
  // BabyTree->Branch("weight_btagsf3_heavy_UP"                                        , &weight_btagsf3_heavy_UP                                                                        );
  // BabyTree->Branch("weight_btagsf3_heavy_DN"                                        , &weight_btagsf3_heavy_DN                                                                        );
  // BabyTree->Branch("weight_btagsf4_heavy_UP"                                        , &weight_btagsf4_heavy_UP                                                                        );
  // BabyTree->Branch("weight_btagsf4_heavy_DN"                                        , &weight_btagsf4_heavy_DN                                                                        );

  BabyTree->Branch("weight_btagsf_iter_central",&weight_btagsf_iter_central);
  BabyTree->Branch("weight_btagsf_iter_down_cferr1",&weight_btagsf_iter_down_cferr1);
  BabyTree->Branch("weight_btagsf_iter_down_cferr2",&weight_btagsf_iter_down_cferr2);
  BabyTree->Branch("weight_btagsf_iter_down_hf",&weight_btagsf_iter_down_hf);
  BabyTree->Branch("weight_btagsf_iter_down_hfstats1",&weight_btagsf_iter_down_hfstats1);
  BabyTree->Branch("weight_btagsf_iter_down_hfstats2",&weight_btagsf_iter_down_hfstats2);
  BabyTree->Branch("weight_btagsf_iter_down_jes",&weight_btagsf_iter_down_jes);
  BabyTree->Branch("weight_btagsf_iter_down_lf",&weight_btagsf_iter_down_lf);
  BabyTree->Branch("weight_btagsf_iter_down_lfstats1",&weight_btagsf_iter_down_lfstats1);
  BabyTree->Branch("weight_btagsf_iter_down_lfstats2",&weight_btagsf_iter_down_lfstats2);
  BabyTree->Branch("weight_btagsf_iter_up_cferr1",&weight_btagsf_iter_up_cferr1);
  BabyTree->Branch("weight_btagsf_iter_up_cferr2",&weight_btagsf_iter_up_cferr2);
  BabyTree->Branch("weight_btagsf_iter_up_hf",&weight_btagsf_iter_up_hf);
  BabyTree->Branch("weight_btagsf_iter_up_hfstats1",&weight_btagsf_iter_up_hfstats1);
  BabyTree->Branch("weight_btagsf_iter_up_hfstats2",&weight_btagsf_iter_up_hfstats2);
  BabyTree->Branch("weight_btagsf_iter_up_jes",&weight_btagsf_iter_up_jes);
  BabyTree->Branch("weight_btagsf_iter_up_lf",&weight_btagsf_iter_up_lf);
  BabyTree->Branch("weight_btagsf_iter_up_lfstats1",&weight_btagsf_iter_up_lfstats1);
  BabyTree->Branch("weight_btagsf_iter_up_lfstats2",&weight_btagsf_iter_up_lfstats2);

  BabyTree->Branch("weight"                                        , &weight                                                                        );
  BabyTree->Branch("weight_lepsf1"                                        , &weight_lepsf1                                                                        );
  BabyTree->Branch("weight_lepsf2"                                        , &weight_lepsf2                                                                        );
  BabyTree->Branch("weight_lepsf3"                                        , &weight_lepsf3                                                                        );
  BabyTree->Branch("weight_lepsf"                                        , &weight_lepsf                                                                        );
  BabyTree->Branch("weight_triggersf"                                        , &weight_triggersf                                                                        );
  BabyTree->Branch("weight_pu"                                        , &weight_pu                                                                        );
  BabyTree->Branch("weight_isrsf"                                        , &weight_isrsf                                                                        );

  BabyTree->Branch("year"                                        , &year                                                                        );
  BabyTree->Branch("yearlumi"                                        , &yearlumi                                                                        );


  // MC syst stuff
  BabyTree->Branch("weight_scale_UP"                                        , &weight_scale_UP                                                                        );
  BabyTree->Branch("weight_scale_DN"                                        , &weight_scale_DN                                                                        );
  BabyTree->Branch("weight_pdf_UP"                                        , &weight_pdf_UP                                                                        );
  BabyTree->Branch("weight_pdf_DN"                                        , &weight_pdf_DN                                                                        );
  BabyTree->Branch("weight_alphas_UP"                                        , &weight_alphas_UP                                                                        );
  BabyTree->Branch("weight_alphas_DN"                                        , &weight_alphas_DN                                                                        );
  BabyTree->Branch("weight_isrvar_UP"                                        , &weight_isrvar_UP                                                                        );
  BabyTree->Branch("weight_isrvar_DN"                                        , &weight_isrvar_DN                                                                        );
  BabyTree->Branch("weight_fsrvar_UP"                                        , &weight_fsrvar_UP                                                                        );
  BabyTree->Branch("weight_fsrvar_DN"                                        , &weight_fsrvar_DN                                                                        );

  // BabyTree->Branch("newweight_scale_UP"  , &newweight_scale_UP );
  // BabyTree->Branch("newweight_scale_DN"  , &newweight_scale_DN );
  // BabyTree->Branch("newweight_pdf_UP"    , &newweight_pdf_UP );
  // BabyTree->Branch("newweight_pdf_DN"    , &newweight_pdf_DN );
  // BabyTree->Branch("newweight_alphas_UP" , &newweight_alphas_UP );
  // BabyTree->Branch("newweight_alphas_DN" , &newweight_alphas_DN );

  //SUSY stuff
  BabyTree->Branch("gl1_p4" , &gl1_p4 );
  BabyTree->Branch("gl2_p4" , &gl2_p4 );
  BabyTree->Branch("glglpt" , &glglpt );
  BabyTree->Branch("isr_unc", &isr_unc);
  BabyTree->Branch("nisrMatch", &nisrMatch);
  BabyTree->Branch("weight_isr"                                           , &weight_isr                                                                           );
  BabyTree->Branch("weight_isr_dy"                                           , &weight_isr_dy                                                                           );
  BabyTree->Branch("weight_isr_tt"                                           , &weight_isr_tt                                                                           );
  BabyTree->Branch("weight_isr_UP"                                        , &weight_isr_UP                                                                        );
  BabyTree->Branch("weight_isr_DN"                                        , &weight_isr_DN                                                                        );

  //Lep3 stuff
  BabyTree->Branch("lep3_el_etaSC"           , &lep3_el_etaSC           );
  BabyTree->Branch("lep3_el_conv_vtx_flag"   , &lep3_el_conv_vtx_flag   );
  BabyTree->Branch("lep3_el_exp_innerlayers" , &lep3_el_exp_innerlayers );
  BabyTree->Branch("lep3_el_threeChargeAgree", &lep3_el_threeChargeAgree);
  BabyTree->Branch("lep3_el_dxyPV"           , &lep3_el_dxyPV           );
  BabyTree->Branch("lep3_el_ip3d"            , &lep3_el_ip3d            );
  BabyTree->Branch("lep3_el_dzPV"            , &lep3_el_dzPV            );
  BabyTree->Branch("lep3_el_MVA_value"       , &lep3_el_MVA_value       );
  BabyTree->Branch("lep3_el_MVA"             , &lep3_el_MVA             );
  // BabyTree->Branch("lep3_iso_RA5"            , &lep3_iso_RA5            );
  BabyTree->Branch("lep3_passes_RA5"         , &lep3_passes_RA5         );
  BabyTree->Branch("lep3_mu_dxyPV"           , &lep3_mu_dxyPV           );
  BabyTree->Branch("lep3_mu_ip3d"            , &lep3_mu_ip3d            );
  BabyTree->Branch("lep3_mu_dzPV"            , &lep3_mu_dzPV            );
  BabyTree->Branch("lep3_mu_ptErr"           , &lep3_mu_ptErr           );

  //Lep 4 stuff
  BabyTree->Branch("lep4_el_etaSC"           , &lep4_el_etaSC           );
  BabyTree->Branch("lep4_el_conv_vtx_flag"   , &lep4_el_conv_vtx_flag   );
  BabyTree->Branch("lep4_el_exp_innerlayers" , &lep4_el_exp_innerlayers );
  BabyTree->Branch("lep4_el_threeChargeAgree", &lep4_el_threeChargeAgree);
  BabyTree->Branch("lep4_el_dxyPV"           , &lep4_el_dxyPV           );
  BabyTree->Branch("lep4_el_ip3d"            , &lep4_el_ip3d            );
  BabyTree->Branch("lep4_el_dzPV"            , &lep4_el_dzPV            );
  BabyTree->Branch("lep4_el_MVA_value"       , &lep4_el_MVA_value       );
  BabyTree->Branch("lep4_el_MVA"             , &lep4_el_MVA             );
  // BabyTree->Branch("lep4_iso_RA5"            , &lep4_iso_RA5            );
  // BabyTree->Branch("lep4_passes_RA5"         , &lep4_passes_RA5         );
  BabyTree->Branch("lep4_mu_dxyPV"           , &lep4_mu_dxyPV           );
  BabyTree->Branch("lep4_mu_ip3d"            , &lep4_mu_ip3d            );
  BabyTree->Branch("lep4_mu_dzPV"            , &lep4_mu_dzPV            );
  BabyTree->Branch("lep4_mu_ptErr"           , &lep4_mu_ptErr           );
  BabyTree->Branch("lep3_isTrigSafeNoIsov1"  , &lep3_isTrigSafeNoIsov1  );
  BabyTree->Branch("lep3_isTrigSafev1"       , &lep3_isTrigSafev1       );
  BabyTree->Branch("lep4_isTrigSafeNoIsov1"  , &lep4_isTrigSafeNoIsov1  );
  BabyTree->Branch("lep4_isTrigSafev1"       , &lep4_isTrigSafev1       );

  BabyTree->Branch("extragenb"       , &extragenb       );
  BabyTree->Branch("extragenbmoms"       , &extragenbmoms       );
  BabyTree->Branch("extragenc"       , &extragenc       );
  BabyTree->Branch("extragencmoms"       , &extragencmoms       );
  BabyTree->Branch("ngenjets"       , &ngenjets       );
  BabyTree->Branch("genht"       , &genht       );
  BabyTree->Branch("ngenjets30"       , &ngenjets30       );
  BabyTree->Branch("genht30"       , &genht30       );
  BabyTree->Branch("passfilter"       , &passfilter       );

  BabyTree->Branch("hasgammatoll"       , &hasgammatoll       );
  BabyTree->Branch("gammatollmomemu"       , &gammatollmomemu       );
  BabyTree->Branch("gammatolldr"       , &gammatolldr       );
  BabyTree->Branch("gammatoll1"       , &gammatoll1       );
  BabyTree->Branch("gammatoll2"       , &gammatoll2       );

  // for tttt
  BabyTree->Branch("bjet_type"       , &bjet_type       );
  BabyTree->Branch("jet_type"       , &jet_type       );
  BabyTree->Branch("ndrlt0p4"       , &ndrlt0p4       );
  BabyTree->Branch("gengood"       , &gengood       );
  BabyTree->Branch("nleptonic"       , &nleptonic       );
  BabyTree->Branch("ntau"       , &ntau       );
  BabyTree->Branch("nleptonicW"       , &nleptonicW       );
  BabyTree->Branch("mfourtop"       , &mfourtop       );
  BabyTree->Branch("nhadronicW"       , &nhadronicW       );
  BabyTree->Branch("nW"       , &nW       );
  BabyTree->Branch("leptonicWSF"       , &leptonicWSF       );
  BabyTree->Branch("hadronicWSF"       , &hadronicWSF       );
  BabyTree->Branch("decayWSF"       , &decayWSF       );
  BabyTree->Branch("njincone1"       , &njincone1       );
  BabyTree->Branch("njincone2"       , &njincone2       );

  // maybe bdt variables?
  BabyTree->Branch("bdt_nforwardjets20", &bdt_nforwardjets20);
  BabyTree->Branch("bdt_avgcdisc", &bdt_avgcdisc);

  // bdt variables
  BabyTree->Branch("bdt_nbtags", &bdt_nbtags);
  BabyTree->Branch("bdt_njets", &bdt_njets);
  BabyTree->Branch("bdt_met", &bdt_met);
  BabyTree->Branch("bdt_ptl2", &bdt_ptl2);
  BabyTree->Branch("bdt_nlb40", &bdt_nlb40);
  BabyTree->Branch("bdt_ntb40", &bdt_ntb40);
  BabyTree->Branch("bdt_nleps", &bdt_nleps);
  BabyTree->Branch("bdt_htb", &bdt_htb);
  BabyTree->Branch("bdt_ml1j1", &bdt_ml1j1);
  BabyTree->Branch("bdt_dphil1l2", &bdt_dphil1l2);
  BabyTree->Branch("bdt_maxmjoverpt", &bdt_maxmjoverpt);
  BabyTree->Branch("bdt_detal1l2", &bdt_detal1l2);
  BabyTree->Branch("bdt_q1", &bdt_q1);
  BabyTree->Branch("bdt_ptj1", &bdt_ptj1);
  BabyTree->Branch("bdt_ptj6", &bdt_ptj6);
  BabyTree->Branch("bdt_ptj7", &bdt_ptj7);
  BabyTree->Branch("bdt_ptj8", &bdt_ptj8);
  BabyTree->Branch("bdt_ptl1", &bdt_ptl1);
  BabyTree->Branch("bdt_ptl3", &bdt_ptl3);
  BabyTree->Branch("bdt_jec_up_nbtags", &bdt_jec_up_nbtags);
  BabyTree->Branch("bdt_jec_dn_nbtags", &bdt_jec_dn_nbtags);
  BabyTree->Branch("bdt_jer_up_nbtags", &bdt_jer_up_nbtags);
  BabyTree->Branch("bdt_jer_dn_nbtags", &bdt_jer_dn_nbtags);
  BabyTree->Branch("bdt_jec_up_njets", &bdt_jec_up_njets);
  BabyTree->Branch("bdt_jec_dn_njets", &bdt_jec_dn_njets);
  BabyTree->Branch("bdt_jer_up_njets", &bdt_jer_up_njets);
  BabyTree->Branch("bdt_jer_dn_njets", &bdt_jer_dn_njets);
  BabyTree->Branch("bdt_jec_up_met", &bdt_jec_up_met);
  BabyTree->Branch("bdt_jec_dn_met", &bdt_jec_dn_met);
  BabyTree->Branch("bdt_jer_up_met", &bdt_jer_up_met);
  BabyTree->Branch("bdt_jer_dn_met", &bdt_jer_dn_met);
  BabyTree->Branch("bdt_jec_up_htb", &bdt_jec_up_htb);
  BabyTree->Branch("bdt_jec_dn_htb", &bdt_jec_dn_htb);
  BabyTree->Branch("bdt_jer_up_htb", &bdt_jer_up_htb);
  BabyTree->Branch("bdt_jer_dn_htb", &bdt_jer_dn_htb);
  BabyTree->Branch("bdt_jec_up_nlb40", &bdt_jec_up_nlb40);
  BabyTree->Branch("bdt_jec_dn_nlb40", &bdt_jec_dn_nlb40);
  BabyTree->Branch("bdt_jer_up_nlb40", &bdt_jer_up_nlb40);
  BabyTree->Branch("bdt_jer_dn_nlb40", &bdt_jer_dn_nlb40);
  BabyTree->Branch("bdt_jec_up_ntb40", &bdt_jec_up_ntb40);
  BabyTree->Branch("bdt_jec_dn_ntb40", &bdt_jec_dn_ntb40);
  BabyTree->Branch("bdt_jer_up_ntb40", &bdt_jer_up_ntb40);
  BabyTree->Branch("bdt_jer_dn_ntb40", &bdt_jer_dn_ntb40);
  BabyTree->Branch("bdt_disc",&bdt_disc);
  BabyTree->Branch("bdt_disc_jec_up",&bdt_disc_jec_up);
  BabyTree->Branch("bdt_disc_jer_up",&bdt_disc_jer_up);
  BabyTree->Branch("bdt_disc_jec_dn",&bdt_disc_jec_dn);
  BabyTree->Branch("bdt_disc_jer_dn",&bdt_disc_jer_dn);
  BabyTree->Branch("bdt_run2_disc",&bdt_run2_disc);
  BabyTree->Branch("bdt_run2_SRDISC",&bdt_run2_SRDISC);
  BabyTree->Branch("bdt_run2_disc_jec_up",&bdt_run2_disc_jec_up);
  BabyTree->Branch("bdt_run2_disc_jer_up",&bdt_run2_disc_jer_up);
  BabyTree->Branch("bdt_run2_disc_jec_dn",&bdt_run2_disc_jec_dn);
  BabyTree->Branch("bdt_run2_disc_jer_dn",&bdt_run2_disc_jer_dn);

  if (applyBtagSFs) {
    // setup btag calibration readers
    // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X
      // calib1 should be period1, calib2 - period2, calib3 - period3, calib4 - INCLUSIVE

    calib1 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_94XSF_V3_B.csv");
    calib2 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_94XSF_V3_C_E.csv");
    calib3 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_94XSF_V3_E_F.csv");
    calib4 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_94XSF_V3_B_F.csv");

    // calib1 = new BTagCalibration("deepflavour", "CORE/Tools/btagsf/data/run2_25ns/DeepFlavour_94XSF_V1_B_F.csv");
    // calib2 = new BTagCalibration("deepflavour", "CORE/Tools/btagsf/data/run2_25ns/DeepFlavour_94XSF_V1_B_F.csv");
    // calib3 = new BTagCalibration("deepflavour", "CORE/Tools/btagsf/data/run2_25ns/DeepFlavour_94XSF_V1_B_F.csv");
    // calib4 = new BTagCalibration("deepflavour", "CORE/Tools/btagsf/data/run2_25ns/DeepFlavour_94XSF_V1_B_F.csv");

    if (gconf.year == 2016) {
        if (gconf.cmssw_ver == 80) {
            calib1 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_Moriond17_B_F.csv"); // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
            calib2 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_Moriond17_G_H.csv"); // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
        } else {
            calib1 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_2016LegacySF_V1.csv"); // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation2016Legacy
        }
    }

    if (gconf.year == 2018) {
        calib1 = new BTagCalibration("deepcsv", "CORE/Tools/btagsf/data/run2_25ns/DeepCSV_102XSF_V1.csv"); // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    }

    // Dylan's bugfix: https://github.com/cmstas/MT2Analysis/commit/04aa9a474da6bccb98abd3a33987f444e6b75e0b#diff-30908e8c1d42c24d6900916966ea4845

    btcr1 = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up","down"});
    btcr1->load(*calib1, BTagEntry::JetFlavor::FLAV_B, "comb");
    btcr1->load(*calib1, BTagEntry::JetFlavor::FLAV_C, "comb");
    btcr1->load(*calib1, BTagEntry::JetFlavor::FLAV_UDSG, "incl");

    btcr2 = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up","down"});
    btcr2->load(*calib2, BTagEntry::JetFlavor::FLAV_B, "comb");
    btcr2->load(*calib2, BTagEntry::JetFlavor::FLAV_C, "comb");
    btcr2->load(*calib2, BTagEntry::JetFlavor::FLAV_UDSG, "incl");

    btcr3 = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up","down"});
    btcr3->load(*calib3, BTagEntry::JetFlavor::FLAV_B, "comb");
    btcr3->load(*calib3, BTagEntry::JetFlavor::FLAV_C, "comb");
    btcr3->load(*calib3, BTagEntry::JetFlavor::FLAV_UDSG, "incl");

    btcr4 = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up","down"});
    btcr4->load(*calib4, BTagEntry::JetFlavor::FLAV_B, "comb");
    btcr4->load(*calib4, BTagEntry::JetFlavor::FLAV_C, "comb");
    btcr4->load(*calib4, BTagEntry::JetFlavor::FLAV_UDSG, "incl");

    btcr_test = new BTagCalibrationReader(BTagEntry::OP_RESHAPING, "central", {"down_cferr1", "down_cferr2", "down_hf", "down_hfstats1", "down_hfstats2", "down_jes", "down_lf", "down_lfstats1", "down_lfstats2", "up_cferr1", "up_cferr2", "up_hf", "up_hfstats1", "up_hfstats2", "up_jes", "up_lf", "up_lfstats1", "up_lfstats2"});
    btcr_test->load(*calib4, BTagEntry::JetFlavor::FLAV_B, "iterativefit");
    btcr_test->load(*calib4, BTagEntry::JetFlavor::FLAV_C, "iterativefit");
    btcr_test->load(*calib4, BTagEntry::JetFlavor::FLAV_UDSG, "iterativefit");

    if (gconf.year == 2016) {
        if (gconf.cmssw_ver == 80) {
            calib_fs = new BTagCalibration("deepcsv_fs", "CORE/Tools/btagsf/data/run2_fastsim/fastsim_deepcsv_ttbar_26_1_2017.csv");
        } else {
            calib_fs = new BTagCalibration("deepcsv_fs", "CORE/Tools/btagsf/data/run2_fastsim/deepcsv_13TEV_16SL_18_3_2019.csv");
        }
    } else if (gconf.year == 2017) {
        calib_fs = new BTagCalibration("deepcsv_fs", "CORE/Tools/btagsf/data/run2_fastsim/deepcsv_13TEV_17SL_18_3_2019.csv");
    } else if (gconf.year == 2018) {
        calib_fs = new BTagCalibration("deepcsv_fs", "CORE/Tools/btagsf/data/run2_fastsim/deepcsv_13TEV_18SL_7_5_2019.csv");
    }

    btcr_fs = new BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up","down"});
    btcr_fs->load(*calib_fs, BTagEntry::JetFlavor::FLAV_B, "fastsim");
    btcr_fs->load(*calib_fs, BTagEntry::JetFlavor::FLAV_C, "fastsim");
    btcr_fs->load(*calib_fs, BTagEntry::JetFlavor::FLAV_UDSG, "fastsim");

    // get btag efficiencies
    TFile* f_btag_eff = 0;
    if (not isFastsim) {

        if (gconf.year == 2016) {
            f_btag_eff = new TFile("CORE/Tools/btagsf/data/run2_25ns/btageff__ttbar_powheg_pythia8_25ns_Moriond17_deepCSV.root");
        } else if (gconf.year == 2017) {
            f_btag_eff = new TFile("CORE/Tools/btagsf/data/run2_25ns/btageff__ttbar_amc_94X_deepCSV.root");
        } else if (gconf.year == 2018) {
            f_btag_eff = new TFile("CORE/Tools/btagsf/data/run2_25ns/btageff__ttbar_amc_102X_deepCSV.root");
        }
        TH2D* h_btag_eff_b_temp    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_b");
        TH2D* h_btag_eff_c_temp    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_c");
        TH2D* h_btag_eff_udsg_temp = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_udsg");

        BabyFile->cd();
        h_btag_eff_b               = (TH2D*) h_btag_eff_b_temp->Clone("h_btag_eff_b");
        h_btag_eff_c               = (TH2D*) h_btag_eff_c_temp->Clone("h_btag_eff_c");
        h_btag_eff_udsg            = (TH2D*) h_btag_eff_udsg_temp->Clone("h_btag_eff_udsg");
        f_btag_eff->Close();
    }
    if (isFastsim) {
        if (gconf.year == 2016) {
            f_btag_eff = new TFile("CORE/Tools/btagsf/data/run2_fastsim/btageff__SMS-T1tttt_2016_80X_deepCSV.root");
        } else if (gconf.year == 2017) {
            f_btag_eff = new TFile("CORE/Tools/btagsf/data/run2_fastsim/btageff__SMS-T1tttt_2017_94X_deepCSV.root");
        } else if (gconf.year == 2018) {
            f_btag_eff = new TFile("CORE/Tools/btagsf/data/run2_fastsim/btageff__DeepCSV_SMS_T2tt_fastsim_Autumn18.root");
        }
        TH2D* h_btag_eff_b_temp    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_b");
        TH2D* h_btag_eff_c_temp    = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_c");
        TH2D* h_btag_eff_udsg_temp = (TH2D*) f_btag_eff->Get("h2_BTaggingEff_csv_med_Eff_udsg");
        BabyFile->cd();
        h_btag_eff_b               = (TH2D*) h_btag_eff_b_temp->Clone("h_btag_eff_b");
        h_btag_eff_c               = (TH2D*) h_btag_eff_c_temp->Clone("h_btag_eff_c");
        h_btag_eff_udsg            = (TH2D*) h_btag_eff_udsg_temp->Clone("h_btag_eff_udsg");
        f_btag_eff->Close();
    }
  }

  // Jet Resolution ... :sad:
  if (gconf.year == 2016) {
      res.loadResolutionFile("CORE/Tools/JetResolution/data/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt");
      res.loadScaleFactorFile("CORE/Tools/JetResolution/data/Summer16_25nsV1_MC_SF_AK4PFchs.txt");
  } else if (gconf.year == 2017) {
      res.loadResolutionFile("CORE/Tools/JetResolution/data/Fall17_V3_MC_PtResolution_AK4PFchs.txt");
      res.loadScaleFactorFile("CORE/Tools/JetResolution/data/Fall17_V3_MC_SF_AK4PFchs.txt");
  } else if (gconf.year == 2018) {
      res.loadResolutionFile("CORE/Tools/JetResolution/data/Fall17_V3_MC_PtResolution_AK4PFchs.txt");
      res.loadScaleFactorFile("CORE/Tools/JetResolution/data/Fall17_V3_MC_SF_AK4PFchs.txt");
  }

  // Load scale1fbs/xsecs from file
  df.loadFromFile("CORE/Tools/datasetinfo/scale1fbs.txt");

}

void babyMaker::InitBabyNtuple(){

    is_fastsim = 0;
    rawmet = -1;
    rawmetPhi = -1;
    calomet = -1;
    calometPhi = -1;
    met = -1;
    metPhi = -1;
    event = -1;
    lumi = -1;
    run = -1;
    filt_csc = 0;
    filt_hbhe = 0;
    filt_hcallaser = 0;
    filt_ecaltp = 0;
    filt_trkfail = 0;
    lep1_isPrompt       = 0;
    lep1_isDirectPrompt = 0;
    lep1_isStat3        = 0;
    lep2_isPrompt       = 0;
    lep2_isDirectPrompt = 0;
    lep2_isStat3        = 0;
    is_real_data = 0;
    scale1fb = 1;
    genps_weight = 1;
    neventstotal = -1;
    xsec_ps = -1;
    kfactor = -1;
    gen_met = -1;
    // genweights.clear();
    pdfweights.clear();
    // genweightsID.clear();
    gen_met_phi = -1;
    skim = 0;
    skim_nomet = 0;
    nleps = 0;
    sr = 0;
    br = 0;
    ss_nleps = -1;
    ss_br = 0;
    ss_sr = -1;
    ss_region = -1;
    njets = -1;
    njetsAG = -1;
    nbtagsAG = -1;
    njets_raw = -1;
    hyp_class = -1;
    ht = -1;
    ht_raw = -1;
    lep1_motherID = 0;
    lep2_motherID = 0;
    lep1_mc_id = -1;
    lep2_mc_id = -1;
    lep1_mc_motherid = -1;
    lep2_mc_motherid = -1;
    lep1_id = -1;
    lep2_id = -1;
    lep1_coneCorrPt = -1;
    lep2_coneCorrPt = -1;
    lep3_coneCorrPt = -1;
    lep4_coneCorrPt = -1;
    lep1_idx = -1;
    lep2_idx = -1;
    jets.clear();
    jets_jec_up.clear();
    jets_jec_dn.clear();
    jets_disc_up.clear();
    jets_disc_dn.clear();
    bjets_jec_up.clear();
    bjets_jec_dn.clear();
    btags_flavor.clear();
    btags_disc.clear();
    btags_cdisc.clear();
    btags_disc_mva.clear();
    btags_disc_ivf.clear();
    ht_raw = -1;
    jets_flavor.clear();
    jets_disc.clear();
    jets_bsf.clear();
    jets_disc_mva.clear();
    jets_disc_ivf.clear();
    jets_JEC.clear();
    btags_JEC.clear();
    jets_undoJEC.clear();
    btags_undoJEC.clear();
    btags.clear();
    jets_unc.clear();
    btags_unc.clear();
    btags_eff.clear();
    btags_effpt.clear();
    btags_sf.clear();
    nbtags = -1;
    nbtags_raw = -1;
    sf_dilepTrig_hpt = -1;
    sf_dilepTrig_lpt = -1;
    sf_dilepTrig_vlpt = -1;
    hyp_type = -1;
    sf_dilep_eff = -1;
    mt = -1;
    mt_l2 = -1;
    mt2 = -1;
    xsec = -1;
    xsec_ps = -1;
    xsec_error = -1;
    mtmin = -1;
    lep3_id = -1;
    lep3_idx = -1;
    lep4_id = -1;
    lep4_idx = -1;
    lep3_quality = -1;
    lep3_mcid = -1;
    lep3_mc_motherid = -1;
    lep3_mc3idx = -1;
    lep3_motherID = -1;
    lep4_mcid = -1;
    lep4_mcidx = -1;
    lep1_iso = -1;
    lep2_iso = -1;
    // lep1_tkIso = -1;
    // lep2_tkIso = -1;
    ncharginos = -1;
    nhiggs = -1;
    higgs_mass = -1;
    genps_p4.clear();
    genjets_p4.clear();
    genps_id.clear();
    genps_id_mother.clear();
    genps_idx_mother.clear();
    genps_status.clear();
    genps_id_grandma.clear();
    all_leps_p4.clear();
    all_leps_id.clear();
    all_leps_veto.clear();
    all_leps_mc3idx.clear();
    lep1_passes_id = false;
    lep2_passes_id = false;
    lep3_passes_id = false;
    lep4_passes_id = false;
    lep1_tight = false;
    lep1_veto = false;
    lep1_fo = false;
    lep2_tight = false;
    lep2_veto = false;
    lep2_fo = false;
    lep3_tight = false;
    lep3_veto = false;
    lep3_fo = false;
    lep4_tight = false;
    lep4_veto = false;
    lep4_fo = false;
    lep1_dxyPV = -999998;
    lep2_dxyPV = -999998;
    lep1_dZ = -999998;
    lep2_dZ = -999998;
    lep1_d0_err = -999998;
    lep2_d0_err = -999998;
    lep1_ip3d = -999998;
    lep2_ip3d = -999998;
    lep1_MVA = -999998;
    lep2_MVA = -999998;
    // lep1_MVA_miniaod = -999998;
    // lep2_MVA_miniaod = -999998;
    // lep1_passes_MVA = false;
    // lep2_passes_MVA = false;
    lep1_ip3d_err = -999998;
    lep2_ip3d_err = -999998;
    lep1_el_conv_vtx_flag = 0;
    lep1_nPixelMiss = 0;
    lep1_el_exp_innerlayers = 0;
    lep1_tightCharge = 0;
    lep1_mu_ptErr = 0;
    lep1_el_threeChargeAgree = 0;
    lep2_el_conv_vtx_flag = 0;
    lep2_nPixelMiss = 0;
    lep2_el_exp_innerlayers = 0;
    lep2_tightCharge = 0;
    lep2_mu_ptErr = 0;
    lep2_el_threeChargeAgree = 0;
    lep1_el_etaSC = -999;
    lep2_el_etaSC = -999;
    lep1_isTrigSafeNoIsov1 = 0;
    lep1_isTrigSafev1 = 0;
    lep2_isTrigSafeNoIsov1 = 0;
    lep2_isTrigSafev1 = 0;
    lep1_tightMuonPOG = 0;
    lep1_mediumMuonPOG = 0;
    lep1_looseMuonPOG = 0;
    lep2_tightMuonPOG = 0;
    lep2_mediumMuonPOG = 0;
    lep2_looseMuonPOG = 0;
    lep3_mediumMuonPOG = 0;

    // FG: add VVV lepton info
    lep1_el_vvv_ss_fo_id=0;
    lep1_el_vvv_3l_fo_id=0;
    lep1_el_vvv_ss_id=0;
    lep1_el_vvv_3l_id=0;
    lep1_mu_vvv_fo_id=0;
    lep1_mu_vvv_id=0;
    lep1_vvv_iso=9999;
    lep2_el_vvv_ss_fo_id=0;
    lep2_el_vvv_3l_fo_id=0;
    lep2_el_vvv_ss_id=0;
    lep2_el_vvv_3l_id=0;
    lep2_mu_vvv_id=0;
    lep2_vvv_iso=9999;

    nVetoElectrons7 = 0;
    nVetoElectrons10 = 0;
    nVetoElectrons25 = 0;
    nVetoMuons5 = 0;
    nVetoMuons10 = 0;
    nVetoMuons25 = 0;
    filename = "";
    lep1_ptrel_ma = -1;
    lep1_ptratio_ma = -1;
    lep2_ptrel_ma = -1;
    lep2_ptratio_ma = -1;
    lep1_ptrel_v1 = -1;
    lep2_ptrel_v1 = -1;
    lep1_ptratio = -1;
    lep2_ptratio = -1;
    lep1_miniIso = -1;
    lep2_miniIso = -1;
    lep1_p4 = LorentzVector(0,0,0,0);
    lep2_p4 = LorentzVector(0,0,0,0);
    lep3_p4 = LorentzVector(0,0,0,0);
    lep4_p4 = LorentzVector(0,0,0,0);
    lep1_p4_gen = LorentzVector(0,0,0,0);
    lep2_p4_gen = LorentzVector(0,0,0,0);
    lep1_closeJet = LorentzVector(0,0,0,0);
    lep2_closeJet = LorentzVector(0,0,0,0);
    dilep_p4 = LorentzVector(0,0,0,0);
    dilep_mass = -1.;
    mass13 = -1.;
    mass23 = -1.;
    mass123 = -1.;
    eleID_kfhits.clear();
    eleID_oldsigmaietaieta.clear();
    eleID_oldsigmaiphiiphi.clear();
    eleID_oldcircularity.clear();
    eleID_oldr9.clear();
    eleID_scletawidth.clear();
    eleID_sclphiwidth.clear();
    eleID_he.clear();
    eleID_psEoverEraw.clear();
    eleID_kfchi2.clear();
    eleID_chi2_hits.clear();
    eleID_fbrem.clear();
    eleID_ep.clear();
    eleID_eelepout.clear();
    eleID_IoEmIop.clear();
    eleID_deltaetain.clear();
    eleID_deltaphiin.clear();
    eleID_deltaetaseed.clear();
    eleID_pT.clear();
    eleID_isbarrel.clear();
    eleID_isendcap.clear();
    eleID_scl_eta.clear();
    muID_dzPV.clear();
    muID_ptSig.clear();
    muID_ip3dSig.clear();
    muID_medMuonPOG.clear();
    mostJets.clear();
    muID_pt.clear();
    muID_eta.clear();
    lep1_isGoodLeg = 0;
    lep2_isGoodLeg = 0;
    lep1_isFakeLeg = 0;
    lep2_isFakeLeg = 0;
    lep1_multiIso          = 0;
    lep2_multiIso          = 0;
    lep1_sip = -1;
    lep2_sip = -1;
    passed_id_inSituFR_lep1 = 0;
    passed_id_inSituFR_lep2 = 0;
    trueNumInt.clear();
    nPUvertices.clear();
    nGoodVertices = 0;
    lep1_trigMatch_noIsoReq = 0;
    lep1_trigMatch_isoReq = 0;
    lep2_trigMatch_noIsoReq = 0;
    lep2_trigMatch_isoReq = 0;
    fired_trigger = 0;
    fired_trigger_ss = 0;
    triggers = 0;
    passes_met_filters = 0;
    passes_chargedcand_filter = 0;
    evt_egclean_pfmet = -1.0;
    evt_muegclean_pfmet = -1.0;
    evt_muegcleanfix_pfmet = -1.0;
    evt_uncorr_pfmet = -1.0;
    filt_noBadMuons = 0;
    filt_duplicateMuons = 0;
    filt_badMuons = 0;
    failsRA2Filter = 0;
    madeExtraZ = 0;
    madeExtraG = 0;
    nisrMatch = -1;
    weight_btagsf  = -9999.;
    weight_btagsf_UP = -9999.;
    weight_btagsf_DN = -9999.;
    weight_btagsf1  = -9999.;
    weight_btagsf1_UP = -9999.;
    weight_btagsf1_DN = -9999.;
    weight_btagsf2  = -9999.;
    weight_btagsf2_UP = -9999.;
    weight_btagsf2_DN = -9999.;
    weight_btagsf3  = -9999.;
    weight_btagsf3_UP = -9999.;
    weight_btagsf3_DN = -9999.;
    weight_btagsf4  = -9999.;
    weight_btagsf4_UP = -9999.;
    weight_btagsf4_DN = -9999.;

    weight_btagsf_light_UP = -9999.;
    weight_btagsf_light_DN = -9999.;
    weight_btagsf1_light_UP = -9999.;
    weight_btagsf1_light_DN = -9999.;
    weight_btagsf2_light_UP = -9999.;
    weight_btagsf2_light_DN = -9999.;
    weight_btagsf3_light_UP = -9999.;
    weight_btagsf3_light_DN = -9999.;
    weight_btagsf4_light_UP = -9999.;
    weight_btagsf4_light_DN = -9999.;
    weight_btagsf_heavy_UP = -9999.;
    weight_btagsf_heavy_DN = -9999.;
    weight_btagsf1_heavy_UP = -9999.;
    weight_btagsf1_heavy_DN = -9999.;
    weight_btagsf2_heavy_UP = -9999.;
    weight_btagsf2_heavy_DN = -9999.;
    weight_btagsf3_heavy_UP = -9999.;
    weight_btagsf3_heavy_DN = -9999.;
    weight_btagsf4_heavy_UP = -9999.;
    weight_btagsf4_heavy_DN = -9999.;

    weight_btagsf_iter_central = -9999.;
    weight_btagsf_iter_down_cferr1 = -9999.;
    weight_btagsf_iter_down_cferr2 = -9999.;
    weight_btagsf_iter_down_hf = -9999.;
    weight_btagsf_iter_down_hfstats1 = -9999.;
    weight_btagsf_iter_down_hfstats2 = -9999.;
    weight_btagsf_iter_down_jes = -9999.;
    weight_btagsf_iter_down_lf = -9999.;
    weight_btagsf_iter_down_lfstats1 = -9999.;
    weight_btagsf_iter_down_lfstats2 = -9999.;
    weight_btagsf_iter_up_cferr1 = -9999.;
    weight_btagsf_iter_up_cferr2 = -9999.;
    weight_btagsf_iter_up_hf = -9999.;
    weight_btagsf_iter_up_hfstats1 = -9999.;
    weight_btagsf_iter_up_hfstats2 = -9999.;
    weight_btagsf_iter_up_jes = -9999.;
    weight_btagsf_iter_up_lf = -9999.;
    weight_btagsf_iter_up_lfstats1 = -9999.;
    weight_btagsf_iter_up_lfstats2 = -9999.;

    weight_isr  = -9999.;
    weight_isr_dy  = -9999.;
    weight_isr_tt  = -9999.;
    weight_isr_UP = -9999.;
    weight_isr_DN = -9999.;
    weight_scale_UP = -9999.;
    weight_scale_DN = -9999.;
    weight_pdf_UP = -9999.;
    weight_pdf_DN = -9999.;
    weight_alphas_UP = -9999.;
    weight_alphas_DN = -9999.;
    // newweight_scale_UP = -9999.;
    // newweight_scale_DN = -9999.;
    // newweight_pdf_UP = -9999.;
    // newweight_pdf_DN = -9999.;
    // newweight_alphas_UP = -9999.;
    // newweight_alphas_DN = -9999.;
    weight_isrvar_UP = -9999.;
    weight_isrvar_DN = -9999.;
    weight_fsrvar_UP = -9999.;
    weight_fsrvar_DN = -9999.;
    mostJets_rawp4.clear();
    mostJets_disc.clear();
    mostJets_unc.clear();
    mostJets_JEC.clear();
    mostJets_undoJEC.clear();
    mostJets_passID.clear();
    mostJets_passCleaning.clear();
    njets_unc_up = 0;
    njets_unc_dn = 0;
    ht_unc_up = 0;
    ht_unc_dn = 0;
    nbtags_unc_up = 0;
    nbtags_unc_dn = 0;
    met_unc_up = 0;
    met_unc_dn = 0;
    metPhi_unc_up = 0;
    metPhi_unc_dn = 0;
    njets_JER_up = 0;
    njets_JER_dn = 0;
    ht_JER_up = 0;
    ht_JER_dn = 0;
    nbtags_JER_up = 0;
    nbtags_JER_dn = 0;
    met_JER_up = 0;
    met_JER_dn = 0;
    metPhi_JER_up = 0;
    metPhi_JER_dn = 0;
    lep1_genps_isHardProcess = 0;
    lep1_genps_fromHardProcessFinalState = 0;
    lep1_genps_fromHardProcessDecayed = 0;
    lep1_genps_isDirectHardProcessTauDecayProductFinalState = 0;
    lep1_genps_fromHardProcessBeforeFSR = 0;
    lep1_genps_isLastCopy = 0;
    lep1_genps_isLastCopyBeforeFSR = 0;
    lep2_genps_isHardProcess = 0;
    lep2_genps_fromHardProcessFinalState = 0;
    lep2_genps_fromHardProcessDecayed = 0;
    lep2_genps_isDirectHardProcessTauDecayProductFinalState = 0;
    lep2_genps_fromHardProcessBeforeFSR = 0;
    lep2_genps_isLastCopy = 0;
    lep2_genps_isLastCopyBeforeFSR = 0;
    glglpt = -1;
    isr_unc = -1;
    gl1_p4 = { 0, 0, 0, 0 };
    gl2_p4 = { 0, 0, 0, 0 };
    lep1_mc3idx = -1;
    lep2_mc3idx = -1;
    lep3_el_etaSC = 0;
    lep3_el_conv_vtx_flag = 0;
    lep3_el_exp_innerlayers = 0;
    lep3_el_threeChargeAgree = 0;
    lep3_el_dxyPV = 0;
    lep3_el_ip3d = 0;
    lep3_el_dzPV = 0;
    lep3_el_MVA_value = 0;
    lep3_el_MVA = 0;
    // lep3_iso_RA5 = 0;
    lep3_passes_RA5 = 0;
    lep3_mu_dxyPV = 0;
    lep3_mu_ip3d = 0;
    lep3_mu_dzPV = 0;
    lep3_mu_ptErr = 0;
    lep4_el_etaSC = 0;
    lep4_el_conv_vtx_flag = 0;
    lep4_el_exp_innerlayers = 0;
    lep4_el_threeChargeAgree = 0;
    lep4_el_dxyPV = 0;
    lep4_el_ip3d = 0;
    lep4_el_dzPV = 0;
    lep4_el_MVA_value = 0;
    lep4_el_MVA = 0;
    // lep4_iso_RA5 = 0;
    // lep4_passes_RA5 = 0;
    lep4_mu_dxyPV = 0;
    lep4_mu_ip3d = 0;
    lep4_mu_dzPV = 0;
    lep4_mu_ptErr = 0;
    lep3_isTrigSafeNoIsov1 = 0;
    lep3_isTrigSafev1 = 0;
    lep4_isTrigSafeNoIsov1 = 0;
    lep4_isTrigSafev1 = 0;
    extragenb = 0;
    extragenbmoms.clear();
    hasgammatoll = false;
    gammatollmomemu = false;
    gammatoll1 = LorentzVector(0,0,0,0);
    gammatoll2 = LorentzVector(0,0,0,0);
    gammatolldr = -1.;
    extragenc = 0;
    extragencmoms.clear();
    ngenjets = 0;
    genht = 0;
    ngenjets30 = 0;
    passfilter = 0;
    genht30 = 0;
    bjet_type.clear();
    jet_type.clear();
    ndrlt0p4 = 0;
    gengood = 0;
    nleptonic = 0;
    ntau = 0;
    nleptonicW = 0;
    mfourtop = 0;
    nhadronicW = 0;
    nW = 0;
    leptonicWSF = 0.;
    hadronicWSF = 0.;
    decayWSF = 0.;
    njincone1 = 0;
    njincone2 = 0;
    prefire_sf = 1.;
    prefire_sfup = 1.;
    prefire_sfdown = 1;
    prefire2016_sf = 1.;
    prefire2016_sfup = 1.;
    prefire2016_sfdown = 1;
    prefire2017_sf = 1.;
    prefire2017_sfup = 1.;
    prefire2017_sfdown = 1;
    prefire2017ele_sf = 1.;
    prefire2017ele_sfup = 1.;
    prefire2017ele_sfdown = 1;
    bdt_nforwardjets20 = 0;
    bdt_avgcdisc = 0.;
    bdt_nbtags = 0.;
    bdt_njets = 0.;
    bdt_met = 0.;
    bdt_ptl2 = 0.;
    bdt_nlb40 = 0.;
    bdt_ntb40 = 0.;
    bdt_nleps = 0.;
    bdt_htb = 0.;
    bdt_ml1j1 = 0.;
    bdt_dphil1l2 = 0.;
    bdt_maxmjoverpt = 0.;
    bdt_detal1l2 = 0.;
    bdt_q1 = 0.;
    bdt_ptj1 = 0.;
    bdt_ptj6 = 0.;
    bdt_ptj7 = 0.;
    bdt_ptj8 = 0.;
    bdt_ptl1 = 0.;
    bdt_ptl3 = 0.;
    bdt_jec_up_nbtags = 0.;
    bdt_jec_dn_nbtags = 0.;
    bdt_jer_up_nbtags = 0.;
    bdt_jer_dn_nbtags = 0.;
    bdt_jec_up_njets = 0.;
    bdt_jec_dn_njets = 0.;
    bdt_jer_up_njets = 0.;
    bdt_jer_dn_njets = 0.;
    bdt_jec_up_met = 0.;
    bdt_jec_dn_met = 0.;
    bdt_jer_up_met = 0.;
    bdt_jer_dn_met = 0.;
    bdt_jec_up_htb = 0.;
    bdt_jec_dn_htb = 0.;
    bdt_jer_up_htb = 0.;
    bdt_jer_dn_htb = 0.;
    bdt_jec_up_nlb40 = 0.;
    bdt_jec_dn_nlb40 = 0.;
    bdt_jer_up_nlb40 = 0.;
    bdt_jer_dn_nlb40 = 0.;
    bdt_jec_up_ntb40 = 0.;
    bdt_jec_dn_ntb40 = 0.;
    bdt_jer_up_ntb40 = 0.;
    bdt_jer_dn_ntb40 = 0.;
    bdt_disc = 0.;
    bdt_disc_jec_up = 0.;
    bdt_disc_jer_up = 0.;
    bdt_disc_jec_dn = 0.;
    bdt_disc_jer_dn = 0.;
    bdt_run2_disc = 0.;
    bdt_run2_SRDISC = -1;
    bdt_run2_disc_jec_up = 0.;
    bdt_run2_disc_jer_up = 0.;
    bdt_run2_disc_jec_dn = 0.;
    bdt_run2_disc_jer_dn = 0.;
    year = -1;
    yearlumi = 0.;
    weight_lepsf1 = 0.;
    weight_lepsf2 = 0.;
    weight_lepsf3 = 0.;
    weight_lepsf = 0.;
    weight_pu = 0.;
    weight_triggersf = 0.;
    weight_isrsf = 0.;
    weight = 0.;
    lep1_pt = -1.;
    lep2_pt = -1.;
    lep3_pt = -1.;
    lep1_eta = -1.;
    lep2_eta = -1.;
    lep3_eta = -1.;
    lep1_phi = -1.;
    lep2_phi = -1.;
    lep3_phi = -1.;
}

//Main function
csErr_t babyMaker::ProcessBaby(string filename_in, FactorizedJetCorrector* jetCorr, JetCorrectionUncertainty* jecUnc, bool isFastsim, int iSignal){


  //Initialize variables
  InitBabyNtuple();

  csErr_t babyErrorStruct;

  //Debug mode
  // if (verbose && evt_cut>0 && tas::evt_event() != evt_cut) return babyErrorStruct;
  // if (evt_cut>0 && tas::evt_event() != evt_cut) return babyErrorStruct;

  //Preliminary stuff
  if (tas::hyp_type().size() < 1) return babyErrorStruct;
  if (tas::mus_dxyPV().size() != tas::mus_dzPV().size()) return babyErrorStruct;

  if (verbose) std::cout << "--------------\nEVENT: " << tas::evt_run() << ":" << tas::evt_lumiBlock() << ":" << tas::evt_event() << "\n--------------" << std::endl;

  //Fill Easy Variables
  filename = filename_in;
  rawmet = evt_pfmet();
  rawmetPhi = evt_pfmetPhi();
  calomet = evt_calomet();
  calometPhi = evt_calometPhi();
  event = tas::evt_event();
  lumi = tas::evt_lumiBlock();
  run = tas::evt_run();
  is_real_data = tas::evt_isRealData();

  // is_miniaodv1 = filename.find("MCRUN2_74_V9") != std::string::npos;
  bool isJetHT = filename.find("JetHT") != std::string::npos;
  bool isRunH = filename.find("Run2016H") != std::string::npos;
  bool isReMiniAOD = filename.find("03Feb2017") != std::string::npos;

  bool get_parentage = true;


  //Fill lepton variables
  hyp_result_t best_hyp_info = chooseBestHyp(verbose);
  hyp_class = best_hyp_info.hyp_class;
  int best_hyp = best_hyp_info.best_hyp;

  // FIXME, don't care about in situ for now, or even nondata OS
  if (hyp_class == 5) return babyErrorStruct;
  if (hyp_class == 7) return babyErrorStruct;
  if (ignore_os) {
      if (!is_real_data && hyp_class == 4) return babyErrorStruct;
  }

  // // if MC and not DY or ttbar, skip opposite sign
  // if ( !is_real_data && !((filename.find("/DY") != std::string::npos) ||
  //             (filename.find("/TT_") != std::string::npos) ||
  //             (filename.find("/TTJets") != std::string::npos) ||
  //             (filename.find("/ST_") != std::string::npos)
  //             ) ) {
  //     if (hyp_class == 4) return babyErrorStruct;
  // }
      

  if (hyp_class == 6){
    madeExtraZ = makesExtraZ(best_hyp).result;
    madeExtraG = makesExtraGammaStar(best_hyp);
  }
  if (verbose) cout << "chose hyp: " << best_hyp << " of class" << hyp_class << endl;
  if (hyp_class == 0 || hyp_class == -1) return babyErrorStruct;


  year = gconf.year;

  if (!is_real_data){
      float sgnMCweight = ((tas::genps_weight() > 0) - (tas::genps_weight() < 0));
      genps_weight = tas::genps_weight();
      trueNumInt = tas::puInfo_trueNumInteractions();
      nPUvertices = puInfo_nPUvertices();
      gen_met = tas::gen_met();
      gen_met_phi = tas::gen_metPhi();

      if (!has_lhe_branches) {

          auto genweights = tas::genweights();  // These two are 20% of the ntuple size
          // genweightsID = tas::genweightsID(); // do we even need the IDs? they are just numbers from 1...N

          //These c-s errors
          int nweights = genweights.size();
          if (nweights>110) {
              float nom = genweights[0];
              // fastsim shifted by 1
              if (isFastsim) nom = genweights[1];
              float sum_pdf = 0.;
              float sum2_pdf = 0.;
              int N = 100;
              for (int ipdf = 9; ipdf < 9+N; ipdf++) {
                  sum_pdf += genweights[ipdf];
                  sum2_pdf += pow(genweights[ipdf],2);
              }
              float rms = sqrt(max(sum2_pdf/N - pow(sum_pdf/N,2),(double)0.0));
              // in 2017, variations are Hessian, so multiply by sqrt(N-1)
              if (gconf.year == 2017 and not isFastsim) {
                  rms *= sqrt(99);
              }
              if (isFastsim) {
                  // fastsim shifted by 1
                  weight_scale_UP = genweights[5] / nom;
                  weight_scale_DN = genweights[9] / nom;
                  weight_alphas_UP = genweights[111] / nom;
                  weight_alphas_DN = genweights[110] / nom;
              } else {
                  weight_scale_UP = genweights[4] / nom;
                  weight_scale_DN = genweights[8] / nom;
                  weight_alphas_UP = genweights[110] / nom;
                  weight_alphas_DN = genweights[109] / nom;
              }
              weight_pdf_UP = (sum_pdf/N+rms) / nom;
              weight_pdf_DN = (sum_pdf/N-rms) / nom;
              weight_fsrvar_UP = 1.;
              weight_fsrvar_DN = 1.;
              weight_isrvar_UP = 1.;
              weight_isrvar_DN = 1.;

              if (fabs(weight_pdf_UP) > 8.) weight_pdf_UP = (weight_pdf_UP > 0 ? 8. : -8.);
              if (fabs(weight_pdf_DN) > 8.) weight_pdf_DN = (weight_pdf_DN > 0 ? 8. : -8.);

              if (read_psweights and nweights > 1000){
                  // 14 extra weights in PSweights samples (2+3*4)
                  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopModGen#Event_Generation
                  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics
                  // go to https://twiki.cern.ch/twiki/pub/CMS/TopModGen/pwg-rwl.txt and grep for <weight id ... there should be ~1080
                  float ps_nom = genweights[nweights-14];
                  if (fabs(ps_nom - nom) < 0.05) {
                      weight_fsrvar_UP = genweights[nweights-11]/ps_nom; // fsr reduced (by sqrt(2) according to recommendation in second link above)
                      weight_fsrvar_DN = genweights[nweights-9]/ps_nom;
                      weight_isrvar_UP = genweights[nweights-8]/ps_nom; // isr default
                      weight_isrvar_DN = genweights[nweights-6]/ps_nom;
                  }
              }

          } else {
              weight_scale_UP = 1.;
              weight_scale_DN = 1.;
              weight_alphas_UP = 1.;
              weight_alphas_DN = 1.;
              weight_pdf_UP = 1.;
              weight_pdf_DN = 1.;
          }

      } else {

          float nom = tas::genHEPMCweight();
          float gen_scale_up = tas::gen_LHEweight_QCDscale_muR2_muF2();
          float gen_scale_dn = tas::gen_LHEweight_QCDscale_muR0p5_muF0p5();
          float gen_pdf_up = tas::gen_LHEweight_PDFVariation_Up();
          float gen_pdf_dn = tas::gen_LHEweight_PDFVariation_Dn();
          float gen_alphas_up = tas::gen_LHEweight_AsMZ_Up();
          float gen_alphas_dn = tas::gen_LHEweight_AsMZ_Dn();

          weight_scale_UP = gen_scale_up;
          weight_scale_DN = gen_scale_dn;
          weight_alphas_UP = gen_alphas_up;
          weight_alphas_DN = gen_alphas_dn;
          weight_pdf_UP = gen_pdf_up;
          weight_pdf_DN = gen_pdf_dn;

          if (fabs(weight_pdf_UP) > 8.) weight_pdf_UP = (weight_pdf_UP > 0 ? 8. : -8.);
          if (fabs(weight_pdf_DN) > 8.) weight_pdf_DN = (weight_pdf_DN > 0 ? 8. : -8.);

          if (read_psweights) {
              auto genweights = tas::genweights();
              int nweights = genweights.size();
              if (nweights > 1000){
                  float ps_nom = genweights[nweights-14];
                  if (fabs(ps_nom - nom) < 0.05) {
                      weight_fsrvar_UP = genweights[nweights-11]/ps_nom; // fsr reduced (by sqrt(2) according to recommendation in second link above)
                      weight_fsrvar_DN = genweights[nweights-9]/ps_nom;
                      weight_isrvar_UP = genweights[nweights-8]/ps_nom; // isr default
                      weight_isrvar_DN = genweights[nweights-6]/ps_nom;
                  } else {
                      weight_fsrvar_UP = 1.;
                      weight_fsrvar_DN = 1.;
                      weight_isrvar_UP = 1.;
                      weight_isrvar_DN = 1.;
                  }
              }
          }



      }
      if (!ignore_scale1fb && iSignal==0) {
          scale1fb = sgnMCweight*df.getScale1fbFromFile(tas::evt_dataset()[0].Data(),tas::evt_CMS3tag()[0].Data());
          xsec = sgnMCweight*df.getXsecFromFile(tas::evt_dataset()[0].Data(),tas::evt_CMS3tag()[0].Data());
          neventstotal = sgnMCweight*df.getnEventsTotalFromFile(tas::evt_dataset()[0].Data(),tas::evt_CMS3tag()[0].Data());
      }
      if (iSignal > 0){
          if (isFastsim) {
              sparms = tas::sparm_values();
          }

          // T1tttt, T5qqqqVV, T5tttt, T5ttcc
          if (iSignal <= 7 or (iSignal==301 and isFastsim)) xsec = go_xsec(sparms[0]).xsec;
          if (iSignal <= 7 or (iSignal==301 and isFastsim)) xsec_error = go_xsec(sparms[0]).percErr;
          if (!isFastsim and ((iSignal == 201) or (iSignal == 301))) {
              // RPV T1qqqqLL, T1tbs
              float mass = -1;
              if (filename.find("mGluino1000") != std::string::npos) mass = 1000;
              else if (filename.find("mGluino1100") != std::string::npos) mass = 1100;
              else if (filename.find("mGluino1200") != std::string::npos) mass = 1200;
              else if (filename.find("mGluino1300") != std::string::npos) mass = 1300;
              else if (filename.find("mGluino1400") != std::string::npos) mass = 1400;
              else if (filename.find("mGluino1500") != std::string::npos) mass = 1500;
              else if (filename.find("mGluino1600") != std::string::npos) mass = 1600;
              else if (filename.find("mGluino1700") != std::string::npos) mass = 1700;
              else if (filename.find("mGluino1800") != std::string::npos) mass = 1800;
              else if (filename.find("mGluino1900") != std::string::npos) mass = 1900;
              else if (filename.find("mGluino2000") != std::string::npos) mass = 2000;
              else if (filename.find("mGluino2100") != std::string::npos) mass = 2100;
              else if (filename.find("mGluino2200") != std::string::npos) mass = 2200;
              else if (filename.find("mGluino2300") != std::string::npos) mass = 2300;
              else if (filename.find("mGluino2400") != std::string::npos) mass = 2400;
              else if (filename.find("mGluino2500") != std::string::npos) mass = 2500;
              else if (filename.find("mGluino2600") != std::string::npos) mass = 2600;
              else {
                  std::cout << "Error, couldn't parse mass from this RPV sample filename" << std::endl;
              }
              xsec = go_xsec(mass).xsec;
              xsec_error = go_xsec(mass).percErr;
              sparms = { mass };
          }
          // T6ttWW, T6ttHZ
          if (iSignal >= 10 and iSignal <= 20) xsec = stop_xsec(sparms[0]).xsec;
          if (iSignal >= 10 and iSignal <= 20) xsec_error = stop_xsec(sparms[0]).percErr;
          // if (iSignal <  100) scale1fb = 1000*xsec/nPoints(iSignal, sparms[0], sparms[1]);
      }
      is_fastsim = isFastsim;
      kfactor = tas::evt_kfactor();

      int higgs_scan = 0;
      if (filename.find("ttH_HToTT") != std::string::npos)  higgs_scan = 1;
      if (filename.find("tHW_HToTT") != std::string::npos)  higgs_scan = 2;
      if (filename.find("tHq_HToTT") != std::string::npos)  higgs_scan = 3;

      if (higgs_scan > 0) {
          for (unsigned int gp=0;gp<tas::genps_id().size();gp++) {
              if (abs(tas::genps_id()[gp])==25) {
                  higgs_mass = floor(tas::genps_p4()[gp].M() + 0.5);
                  break;
              }
          }
          xsec = xsec_higgs(higgs_scan, higgs_mass);
          xsec_ps = xsec_higgs(higgs_scan+3, higgs_mass);
          if (iSignal > 0) {
              scale1fb = 1000*xsec/nPoints_higgs(higgs_scan, higgs_mass);
          }
      }

  }

  lep1_p4 = (tas::hyp_ll_p4().at(best_hyp).pt() > tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_p4().at(best_hyp) : tas::hyp_lt_p4().at(best_hyp);
  lep2_p4 = (tas::hyp_ll_p4().at(best_hyp).pt() <= tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_p4().at(best_hyp) : tas::hyp_lt_p4().at(best_hyp);
  lep1_id = (tas::hyp_ll_p4().at(best_hyp).pt() > tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_id().at(best_hyp) : tas::hyp_lt_id().at(best_hyp);
  lep2_id = (tas::hyp_ll_p4().at(best_hyp).pt() <= tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_id().at(best_hyp) : tas::hyp_lt_id().at(best_hyp);
  lep1_idx = (tas::hyp_ll_p4().at(best_hyp).pt() > tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_index().at(best_hyp) : tas::hyp_lt_index().at(best_hyp);
  lep2_idx = (tas::hyp_ll_p4().at(best_hyp).pt() <= tas::hyp_lt_p4().at(best_hyp).pt()) ? tas::hyp_ll_index().at(best_hyp) : tas::hyp_lt_index().at(best_hyp);
  lep1_coneCorrPt = coneCorrPt(lep1_id, lep1_idx);
  lep2_coneCorrPt = coneCorrPt(lep2_id, lep2_idx);
  Lep lep1 = Lep(lep1_id, lep1_idx);
  Lep lep2 = Lep(lep2_id, lep2_idx);
  lep1_dxyPV = lep1.dxyPV();
  lep2_dxyPV = lep2.dxyPV();
  lep1_dZ = lep1.dzPV();
  lep2_dZ = lep2.dzPV();
  lep1_d0_err = lep1.d0Err();
  lep2_d0_err = lep2.d0Err();
  lep1_ip3d = lep1.ip3d();
  lep2_ip3d = lep2.ip3d();
  lep1_ip3d_err = lep1.ip3dErr();
  lep2_ip3d_err = lep2.ip3dErr();
  hyp_type = tas::hyp_type().at(best_hyp);
  pair <Lep, int> thirdLepton = make_pair(Lep(0,-1), -1);
  bool dolowpt = true;
  if (hyp_class == 8) {
      // thirdLepton =  make_pair(best_lep3,2);
  } else {
      thirdLepton = getThirdLepton(best_hyp, dolowpt);
  }
  lep3_id = thirdLepton.first.pdgId();
  lep3_idx = thirdLepton.first.idx();
  if (lep3_idx >= 0 && (abs(lep3_id) == 11 || abs(lep3_id) == 13)){
    lep3_p4 = thirdLepton.first.p4();
    lep3_coneCorrPt = coneCorrPt(lep3_id, lep3_idx);
    if (!is_real_data && get_parentage){
      pair <int, int> lep3_parentage = lepMotherID_v2(thirdLepton.first);
      lep3_motherID = lep3_parentage.first;
      lep3_mc3idx = lep3_parentage.second;
      lep3_mcid = lep3_mc3idx >= 0 ? tas::genps_id().at(lep3_mc3idx) : -9999;
      lep3_mc_motherid                                              = (lep3_mc3idx < 0) ? 0 : tas::genps_id_mother().at(lep3_mc3idx);
    }
  }
  lep3_quality = thirdLepton.second;
  Lep fourthLepton = getFourthLepton(best_hyp, lep3_id, lep3_idx);
  // // same function for getting 3rd lepton, but now we skip over the 3rd lepton we already found
  // Lep fourthLepton = getThirdLepton(best_hyp, lep3_id, lep3_idx).first;
  lep4_id = fourthLepton.pdgId();
  lep4_idx = fourthLepton.idx();
  if (lep4_idx >= 0 && (abs(lep4_id) == 11 || abs(lep4_id) == 13)){
    lep4_p4 = fourthLepton.p4();
    lep4_coneCorrPt = coneCorrPt(lep4_id, lep4_idx);
    if (!is_real_data && get_parentage){
      pair <int, int> lep4_parentage = lepMotherID_v2(fourthLepton);
      lep4_mcidx = lep4_parentage.second;
      lep4_mcid = lep4_mcidx >= 0 ? tas::genps_id().at(lep4_mcidx) : -9999;
    }
  }
  lep1_iso = abs(lep1_id) == 11 ? eleRelIso03(lep1_idx, SS) :  muRelIso03(lep1_idx, SS);
  lep2_iso = abs(lep2_id) == 11 ? eleRelIso03(lep2_idx, SS) :  muRelIso03(lep2_idx, SS);
  // lep1_tkIso = abs(lep1_id) == 11 ? els_tkIso().at(lep1_idx)/lep1_p4.pt() : mus_iso03_sumPt().at(lep1_idx)/lep1_p4.pt();
  // lep2_tkIso = abs(lep2_id) == 11 ? els_tkIso().at(lep2_idx)/lep2_p4.pt() : mus_iso03_sumPt().at(lep2_idx)/lep2_p4.pt();
  lep1_multiIso = abs(lep1_id) == 11 ? passMultiIso(11, lep1_idx, gconf.multiiso_el_minireliso, gconf.multiiso_el_ptratio, gconf.multiiso_el_ptrel, gconf.ea_version, 2) : passMultiIso(13, lep1_idx, gconf.multiiso_mu_minireliso, gconf.multiiso_mu_ptratio, gconf.multiiso_mu_ptrel, gconf.ea_version, 2);
  lep2_multiIso = abs(lep2_id) == 11 ? passMultiIso(11, lep2_idx, gconf.multiiso_el_minireliso, gconf.multiiso_el_ptratio, gconf.multiiso_el_ptrel, gconf.ea_version, 2) : passMultiIso(13, lep2_idx, gconf.multiiso_mu_minireliso, gconf.multiiso_mu_ptratio, gconf.multiiso_mu_ptrel, gconf.ea_version, 2);
  lep1_sip = abs(lep1_id) == 11 ? fabs(els_ip3d().at(lep1_idx))/els_ip3derr().at(lep1_idx) : fabs(mus_ip3d().at(lep1_idx))/mus_ip3derr().at(lep1_idx);
  lep2_sip = abs(lep2_id) == 11 ? fabs(els_ip3d().at(lep2_idx))/els_ip3derr().at(lep2_idx) : fabs(mus_ip3d().at(lep2_idx))/mus_ip3derr().at(lep2_idx);
  dilep_p4 = lep1_p4 + lep2_p4;
  dilep_mass = dilep_p4.M();
  lep1_passes_id = isGoodLepton(lep1_id, lep1_idx);
  lep2_passes_id = isGoodLepton(lep2_id, lep2_idx);
  lep1_MVA = abs(lep1_id) == 11 ? getMVAoutput(lep1_idx,true) : -9999;
  lep2_MVA = abs(lep2_id) == 11 ? getMVAoutput(lep2_idx,true) : -9999;
  // lep1_MVA_miniaod = lep1_MVA;
  // lep2_MVA_miniaod = lep2_MVA;

  // lep1_passes_MVA = abs(lep1_id) == 11 ? passesMVAforID(lep1_idx, SS_medium_noiso_v6) : true;
  // lep2_passes_MVA = abs(lep2_id) == 11 ? passesMVAforID(lep2_idx, SS_medium_noiso_v6) : true;

  // Reject events that fail trigger matching
  bool isLLcc = lep1_coneCorrPt < 25. and lep2_coneCorrPt < 25.;
  if (gconf.year == 2016) {
      if (ht < 300 && hyp_type != 0){
          if (abs(lep1_id) == 11 && !isTriggerSafe_v1(lep1_idx)) return babyErrorStruct;
          if (abs(lep2_id) == 11 && !isTriggerSafe_v1(lep2_idx)) return babyErrorStruct;
      }
  } else if (gconf.year == 2017) {
      if (hyp_type != 0 and !isLLcc){
          if (abs(lep1_id) == 11 && !isTriggerSafe_v1(lep1_idx)) return babyErrorStruct;
          if (abs(lep2_id) == 11 && !isTriggerSafe_v1(lep2_idx)) return babyErrorStruct;
      }
  } else if (gconf.year == 2018) {
      if (hyp_type != 0 and !isLLcc){
          if (abs(lep1_id) == 11 && !isTriggerSafe_v1(lep1_idx)) return babyErrorStruct;
          if (abs(lep2_id) == 11 && !isTriggerSafe_v1(lep2_idx)) return babyErrorStruct;
      }
  }
  if (verbose) cout << "passed trigger safety" << endl;

  //Corrected MET
  int use_cleaned_met = 0;
  if (tas::evt_isRealData() and (filename.find("-03Feb2017-") != std::string::npos)) {
      use_cleaned_met = 1;
  }
  if (gconf.year == 2017) {
      use_cleaned_met = 2;
  }
  pair <float, float> T1CHSMET = getT1CHSMET_fromMINIAOD(jetCorr, NULL, 0, false, use_cleaned_met);
  met = T1CHSMET.first;
  metPhi = T1CHSMET.second;

  //More First lepton stuff
  if (abs(lep1_id) == 11 || abs(lep1_id) == 13){
    lep1_passes_id = isGoodLepton(lep1_id, lep1_idx);
    lep1_tight = abs(lep1_id) == 11 ? isGoodElectron(lep1_idx) : isGoodMuon(lep1_idx);
    lep1_veto = abs(lep1_id) == 11 ? isGoodVetoElectron(lep1_idx) : isGoodVetoMuon(lep1_idx);
    lep1_fo = abs(lep1_id) == 11 ? isFakableElectron(lep1_idx) : isFakableMuon(lep1_idx);
    lep1_el_conv_vtx_flag = abs(lep1_id) == 11 ? (!(els_conv_vtx_flag().at(lep1_idx))) : 0;
    lep1_nPixelMiss = abs(lep1_id) == 11 ? (!(els_lost_pixelhits().at(lep1_idx))) : -1;
    lep1_tightCharge = abs(lep1_id) == 11 ? tightChargeEle(lep1_idx) : -1;
    lep1_el_exp_innerlayers  = abs(lep1_id) == 11 ? (els_exp_innerlayers().at(lep1_idx)) : -1;
    lep1_el_threeChargeAgree = abs(lep1_id) == 11 ? threeChargeAgree(lep1_idx) : -1;
    lep1_el_etaSC = abs(lep1_id) == 11 ? els_etaSC().at(lep1_idx) : -999;

    lep1_mu_ptErr = abs(lep1_id) == 13 ? (mus_ptErr().at(lep1_idx)/mus_trk_p4().at(lep1_idx).pt() < 0.2) : -1;

    lep1_isTrigSafeNoIsov1   = abs(lep1_id) == 11 ? isTriggerSafenoIso_v1(lep1_idx) : 1;
    lep1_isTrigSafev1        = abs(lep1_id) == 11 ? isTriggerSafe_v1(lep1_idx) : 1;

    lep1_tightMuonPOG = abs(lep1_id) == 13 ? isTightMuonPOG(lep1_idx) : 0;
    lep1_mediumMuonPOG = abs(lep1_id) == 13 ? isMediumMuonPOG(lep1_idx) : 0;
    lep1_looseMuonPOG = abs(lep1_id) == 13 ? isLooseMuonPOG(lep1_idx) : 0;

    lep1_el_vvv_ss_fo_id = abs(lep1_id) == 11 ? electronID(lep1_idx,VVV_fo_noiso_v5) : 0;
    lep1_el_vvv_3l_fo_id = abs(lep1_id) == 11 ? electronID(lep1_idx,VVV_3l_fo_noiso_v5) : 0;
    lep1_el_vvv_ss_id = abs(lep1_id) == 11 ? electronID(lep1_idx,VVV_tight_noiso_v5) : 0;
    lep1_el_vvv_3l_id = abs(lep1_id) == 11 ? electronID(lep1_idx,VVV_3l_tight_noiso_v5) : 0;
    lep1_mu_vvv_fo_id = abs(lep1_id) == 13 ? muonID(lep1_idx,VVV_fo_noiso_v5) : 0;
    lep1_mu_vvv_id = abs(lep1_id) == 13 ? muonID(lep1_idx,VVV_tight_noiso_v5) : 0;
    lep1_vvv_iso = abs(lep1_id) == 13 ? (gconf.year == 2016 ? muRelIso03EA(lep1_idx, 2, true) : muRelIso03EA(lep1_idx,4,true)) : abs(lep1_id) == 11 ? (gconf.year == 2016 ? eleRelIso03EA(lep1_idx,2,true) : eleRelIso03EA(lep1_idx,4,true)) : 0;
  }


  //More Second lepton stuff
  if (abs(lep2_id) == 11 || abs(lep2_id) == 13){
    lep2_passes_id = isGoodLepton(lep2_id, lep2_idx);
    lep2_tight = abs(lep2_id) == 11 ? isGoodElectron(lep2_idx) : isGoodMuon(lep2_idx);
    lep2_veto = abs(lep2_id) == 11 ? isGoodVetoElectron(lep2_idx) : isGoodVetoMuon(lep2_idx);
    lep2_fo = abs(lep2_id) == 11 ? isFakableElectron(lep2_idx) : isFakableMuon(lep2_idx);
    lep2_el_conv_vtx_flag = abs(lep2_id) == 11 ? (!(els_conv_vtx_flag().at(lep2_idx))) : 0;
    lep2_nPixelMiss = abs(lep2_id) == 11 ? (!(els_lost_pixelhits().at(lep2_idx))) : -1;
    lep2_tightCharge = abs(lep2_id) == 11 ? tightChargeEle(lep2_idx) : -1;
    lep2_el_exp_innerlayers  = abs(lep2_id) == 11 ? (els_exp_innerlayers().at(lep2_idx)) : -1;
    lep2_el_threeChargeAgree = abs(lep2_id) == 11 ? threeChargeAgree(lep2_idx) : -1;
    lep2_el_etaSC = abs(lep2_id) == 11 ? els_etaSC().at(lep2_idx) : -999;

    lep2_mu_ptErr = abs(lep2_id) == 13 ? (mus_ptErr().at(lep2_idx)/mus_trk_p4().at(lep2_idx).pt() < 0.2) : -1;

    lep2_isTrigSafeNoIsov1   = abs(lep2_id) == 11 ? isTriggerSafenoIso_v1(lep2_idx) : 1;
    lep2_isTrigSafev1        = abs(lep2_id) == 11 ? isTriggerSafe_v1(lep2_idx) : 1;

    lep2_tightMuonPOG = abs(lep2_id) == 13 ? isTightMuonPOG(lep2_idx) : 0;
    lep2_mediumMuonPOG = abs(lep2_id) == 13 ? isMediumMuonPOG(lep2_idx) : 0;
    lep2_looseMuonPOG = abs(lep2_id) == 13 ? isLooseMuonPOG(lep2_idx) : 0;

    lep2_el_vvv_ss_fo_id = abs(lep2_id) == 11 ? electronID(lep2_idx,VVV_fo_noiso_v5) : 0;
    lep2_el_vvv_3l_fo_id = abs(lep2_id) == 11 ? electronID(lep2_idx,VVV_3l_fo_noiso_v5) : 0;
    lep2_el_vvv_ss_id = abs(lep2_id) == 11 ? electronID(lep2_idx,VVV_tight_noiso_v5) : 0;
    lep2_el_vvv_3l_id = abs(lep2_id) == 11 ? electronID(lep2_idx,VVV_3l_tight_noiso_v5) : 0;
    lep2_mu_vvv_fo_id = abs(lep2_id) == 13 ? muonID(lep2_idx,VVV_fo_noiso_v5) : 0;
    lep2_mu_vvv_id = abs(lep2_id) == 13 ? muonID(lep2_idx,VVV_tight_noiso_v5) : 0;
    lep2_vvv_iso = abs(lep2_id) == 13 ? (gconf.year == 2016 ? muRelIso03EA(lep2_idx, 2, true) : muRelIso03EA(lep2_idx,4,true)) : abs(lep2_id) == 11 ? (gconf.year == 2016 ? eleRelIso03EA(lep2_idx,2,true) : eleRelIso03EA(lep2_idx,4,true)) : 0;
  }

  //More Third lepton stuff
  if (abs(lep3_id) == 11 || abs(lep3_id) == 13){
    lep3_passes_id = isGoodLepton(lep3_id, lep3_idx);
    lep3_tight = abs(lep3_id) == 11 ? isGoodElectron(lep3_idx) : isGoodMuon(lep3_idx);
    lep3_veto = abs(lep3_id) == 11 ? isGoodVetoElectron(lep3_idx) : isGoodVetoMuon(lep3_idx);
    lep3_fo = abs(lep3_id) == 11 ? isFakableElectron(lep3_idx) : isFakableMuon(lep3_idx);
    if (lep3_veto || hyp_class==8){
      if (abs(lep3_id) == 11){

          float ptmin = 15;
          float ptmax = 25;
          float elpt = els_p4().at(lep3_idx).pt();
          auto mvacut = [ptmin,ptmax](float A, float B, float C, float pt_) {
              return pt_>10 ? std::min(A, std::max(B,A + (B-A)/(ptmax-ptmin)*(pt_-ptmin))) : C;
          };

        float etaSC              = fabs(els_etaSC().at(lep3_idx));
        lep3_el_etaSC            = (etaSC <= 2.5);
        lep3_el_conv_vtx_flag    = (!(els_conv_vtx_flag().at(lep3_idx)));
        lep3_el_exp_innerlayers  = (els_exp_innerlayers().at(lep3_idx) <= 0);
        lep3_el_threeChargeAgree = (threeChargeAgree(lep3_idx));
        lep3_el_dxyPV            = (fabs(els_dxyPV().at(lep3_idx)) <= 0.05);
        lep3_el_ip3d             = (fabs(els_ip3d().at(lep3_idx))/els_ip3derr().at(lep3_idx) < 4);
        lep3_el_dzPV             = (fabs(els_dzPV().at(lep3_idx)) < 0.1);
        lep3_el_MVA_value        = getMVAoutput(lep3_idx,true);
        lep3_el_MVA = (etaSC < 0.8) ? (lep3_el_MVA_value > mvacut(0.77,0.52,0.77,elpt)) : ((etaSC >= 0.8 && etaSC <= 1.479) ? lep3_el_MVA_value > mvacut(0.56,0.11,0.56,elpt) : lep3_el_MVA_value > mvacut(0.48,-0.01,0.48,elpt));
        // lep3_iso_RA5             = passMultiIso(11, lep3_idx, 0.12, 0.80, 7.2, gconf.ea_version, 2);
        lep3_isTrigSafeNoIsov1   = isTriggerSafenoIso_v1(lep3_idx);
        lep3_isTrigSafev1        = isTriggerSafe_v1(lep3_idx);
        lep3_passes_RA5 = (lep3_el_etaSC && lep3_el_conv_vtx_flag && lep3_el_exp_innerlayers && lep3_el_threeChargeAgree && lep3_el_dxyPV && lep3_el_ip3d && lep3_el_dzPV && lep3_el_MVA && lep3_iso_RA5 && lep3_isTrigSafeNoIsov1);
      }
      if (abs(lep3_id) == 13){
        lep3_mu_dxyPV = (fabs(mus_dxyPV().at(lep3_idx)) <= 0.05);
        lep3_mu_ip3d  = (fabs(mus_ip3d().at(lep3_idx))/mus_ip3derr().at(lep3_idx) < 4);
        lep3_mu_dzPV  = (fabs(mus_dzPV().at(lep3_idx)) <= 0.1);
        lep3_mu_ptErr = (mus_ptErr().at(lep3_idx)/mus_trk_p4().at(lep3_idx).pt() < 0.2);
        // lep3_iso_RA5  = passMultiIso(13, lep3_idx, 0.16, 0.76, 7.2, gconf.ea_version, 2);
        lep3_mediumMuonPOG = isMediumMuonPOG(lep3_idx);
        lep3_passes_RA5 = (lep3_mu_dxyPV && lep3_mu_ip3d && lep3_mu_dzPV && lep3_mu_ptErr && lep3_iso_RA5 && lep3_mediumMuonPOG);
      }
    }
  }

  //More fourth lepton stuff
  if (abs(lep4_id) == 11 || abs(lep4_id) == 13){
    lep4_veto = abs(lep4_id) == 11 ? isGoodVetoElectron(lep4_idx) : isGoodVetoMuon(lep4_idx);
    if (lep4_veto){
      lep4_passes_id = isGoodLepton(lep4_id, lep4_idx);
      lep4_tight = abs(lep4_id) == 11 ? isGoodElectron(lep4_idx) : isGoodMuon(lep4_idx);
      lep4_fo = abs(lep4_id) == 11 ? isFakableElectron(lep4_idx) : isFakableMuon(lep4_idx);
      if (abs(lep4_id) == 11){
        float etaSC              = fabs(els_etaSC().at(lep4_idx));
        lep4_el_etaSC            = (etaSC <= 2.5);
        lep4_el_conv_vtx_flag    = (!(els_conv_vtx_flag().at(lep4_idx)));
        lep4_el_exp_innerlayers  = (els_exp_innerlayers().at(lep4_idx) <= 0);
        lep4_el_threeChargeAgree = (threeChargeAgree(lep4_idx));
        lep4_el_dxyPV            = (fabs(els_dxyPV().at(lep4_idx)) <= 0.05);
        lep4_el_ip3d             = (fabs(els_ip3d().at(lep4_idx))/els_ip3derr().at(lep4_idx) < 4);
        lep4_el_dzPV             = (fabs(els_dzPV().at(lep4_idx)) < 0.1);
        lep4_el_MVA_value        = getMVAoutput(lep4_idx,true);
        lep4_el_MVA              = ((etaSC < 0.8) ? (lep4_el_MVA_value > 0.87) : ((etaSC <= 1.479) ? (lep4_el_MVA_value > 0.60) : (lep4_el_MVA_value > 0.17)));
        // lep4_iso_RA5             = passMultiIso(11, lep4_idx, 0.16, 0.76, 7.2, gconf.ea_version, 2);
        lep4_isTrigSafeNoIsov1   = isTriggerSafenoIso_v1(lep4_idx);
        lep4_isTrigSafev1        = isTriggerSafe_v1(lep4_idx);
        // lep4_passes_RA5 = (lep4_el_etaSC && lep4_el_conv_vtx_flag && lep4_el_exp_innerlayers && lep4_el_threeChargeAgree && lep4_el_dxyPV && lep4_el_ip3d && lep4_el_dzPV && lep4_el_MVA && lep4_iso_RA5 && lep4_isTrigSafeNoIsov1 && (ht > 400 || lep4_isTrigSafev1));
      }
      if (abs(lep4_id) == 13){
        lep4_mu_dxyPV = (fabs(mus_dxyPV().at(lep4_idx)) <= 0.05);
        lep4_mu_ip3d  = (fabs(mus_ip3d().at(lep4_idx))/mus_ip3derr().at(lep4_idx) < 4);
        lep4_mu_dzPV  = (fabs(mus_dzPV().at(lep4_idx)) <= 0.1);
        lep4_mu_ptErr = (mus_ptErr().at(lep4_idx)/mus_trk_p4().at(lep4_idx).pt() < 0.2);
        // lep4_iso_RA5  = passMultiIso(13, lep4_idx, 0.12, 0.80, 7.2, gconf.ea_version, 2);
        // lep4_passes_RA5 = (lep4_mu_dxyPV && lep4_mu_ip3d && lep4_mu_dzPV && lep4_mu_ptErr && lep4_iso_RA5);
      }
    }
  }

  //Lepton MC variables
  if (!is_real_data && get_parentage){
    pair <int, int> lep1_parentage                          = lepMotherID_v2(lep1);
    pair <int, int> lep2_parentage                          = lepMotherID_v2(lep2);
    lep1_mc3idx                                             = lep1_parentage.second;
    lep2_mc3idx                                             = lep2_parentage.second;
    lep1_mc_id                                              = (lep1_mc3idx < 0) ? 0 : tas::genps_id().at(lep1_mc3idx);
    lep2_mc_id                                              = (lep2_mc3idx < 0) ? 0 : tas::genps_id().at(lep2_mc3idx);
    lep1_p4_gen                                              = (lep1_mc3idx < 0) ? LorentzVector(0,0,0,0) : tas::genps_p4().at(lep1_mc3idx);
    lep2_p4_gen                                              = (lep2_mc3idx < 0) ? LorentzVector(0,0,0,0) : tas::genps_p4().at(lep2_mc3idx);
    lep1_mc_motherid                                              = (lep1_mc3idx < 0) ? 0 : tas::genps_id_mother().at(lep1_mc3idx);
    lep2_mc_motherid                                              = (lep2_mc3idx < 0) ? 0 : tas::genps_id_mother().at(lep2_mc3idx);
    lep1_isPrompt                                           = (lep1_mc3idx < 0) ? 0 : tas::genps_isPromptFinalState().at(lep1_mc3idx);
    // lep1_isStat3                                            = (lep1_mc3idx < 0) ? 0 : tas::genps_isMostlyLikePythia6Status3().at(lep1_mc3idx);
    lep1_isDirectPrompt                                     = (lep1_mc3idx < 0) ? 0 : tas::genps_isDirectPromptTauDecayProductFinalState().at(lep1_mc3idx);
    lep1_genps_isHardProcess                                = (lep1_mc3idx < 0) ? 0 : tas::genps_isHardProcess().at(lep1_mc3idx);
    lep1_genps_fromHardProcessFinalState                    = (lep1_mc3idx < 0) ? 0 : tas::genps_fromHardProcessFinalState().at(lep1_mc3idx);
    lep1_genps_fromHardProcessDecayed                       = (lep1_mc3idx < 0) ? 0 : tas::genps_fromHardProcessDecayed().at(lep1_mc3idx);
    lep1_genps_isDirectHardProcessTauDecayProductFinalState = (lep1_mc3idx < 0) ? 0 : tas::genps_isDirectHardProcessTauDecayProductFinalState().at(lep1_mc3idx);
    lep1_genps_fromHardProcessBeforeFSR                     = (lep1_mc3idx < 0) ? 0 : tas::genps_fromHardProcessBeforeFSR().at(lep1_mc3idx);
    lep1_genps_isLastCopy                                   = (lep1_mc3idx < 0) ? 0 : tas::genps_isLastCopy().at(lep1_mc3idx);
    lep1_genps_isLastCopyBeforeFSR                          = (lep1_mc3idx < 0) ? 0 : tas::genps_isLastCopyBeforeFSR().at(lep1_mc3idx);
    lep1_motherID                                           = lep1_parentage.first;
    lep2_motherID                                           = lep2_parentage.first;
    lep2_isPrompt                                           = (lep2_mc3idx < 0) ? 0 : tas::genps_isPromptFinalState().at(lep2_mc3idx);
    lep2_isDirectPrompt                                     = (lep2_mc3idx < 0) ? 0 : tas::genps_isDirectPromptTauDecayProductFinalState().at(lep2_mc3idx);
    // lep2_isStat3                                            = (lep2_mc3idx < 0) ? 0 : tas::genps_isMostlyLikePythia6Status3().at(lep2_mc3idx);
    lep2_genps_isHardProcess                                = (lep2_mc3idx < 0) ? 0 : tas::genps_isHardProcess().at(lep2_mc3idx);
    lep2_genps_fromHardProcessFinalState                    = (lep2_mc3idx < 0) ? 0 : tas::genps_fromHardProcessFinalState().at(lep2_mc3idx);
    lep2_genps_fromHardProcessDecayed                       = (lep2_mc3idx < 0) ? 0 : tas::genps_fromHardProcessDecayed().at(lep2_mc3idx);
    lep2_genps_isDirectHardProcessTauDecayProductFinalState = (lep2_mc3idx < 0) ? 0 : tas::genps_isDirectHardProcessTauDecayProductFinalState().at(lep2_mc3idx);
    lep2_genps_fromHardProcessBeforeFSR                     = (lep2_mc3idx < 0) ? 0 : tas::genps_fromHardProcessBeforeFSR().at(lep2_mc3idx);
    lep2_genps_isLastCopy                                   = (lep2_mc3idx < 0) ? 0 : tas::genps_isLastCopy().at(lep2_mc3idx);
    lep2_genps_isLastCopyBeforeFSR                          = (lep2_mc3idx < 0) ? 0 : tas::genps_isLastCopyBeforeFSR().at(lep2_mc3idx);
  }


  //PtRel for both leptons
  lep1_ptrel_v1 = getPtRel(lep1_id, lep1_idx, true, ssWhichCorr);
  lep2_ptrel_v1 = getPtRel(lep2_id, lep2_idx, true, ssWhichCorr);

  //Closest jet for both leptons
  lep1_closeJet = closestJet(lep1_p4, 0.4, 3.0, ssWhichCorr);
  lep2_closeJet = closestJet(lep2_p4, 0.4, 3.0, ssWhichCorr);

  lep1_ptratio = lep1_p4.pt()/lep1_closeJet.pt();
  lep2_ptratio = lep2_p4.pt()/lep2_closeJet.pt();

  //MiniIso
  lep1_miniIso = abs(lep1_id)==11 ? elMiniRelIsoCMS3_EA(lep1_idx, gconf.ea_version) : muMiniRelIsoCMS3_EA(lep1_idx, gconf.ea_version);
  lep2_miniIso = abs(lep2_id)==11 ? elMiniRelIsoCMS3_EA(lep2_idx, gconf.ea_version) : muMiniRelIsoCMS3_EA(lep2_idx, gconf.ea_version);

  //For inSituFR, both must pass looser ID (easier than selection ID)
  passed_id_inSituFR_lep1 = isInSituFRLepton(lep1_id, lep1_idx);
  passed_id_inSituFR_lep2 = isInSituFRLepton(lep2_id, lep2_idx);



  // This is 50% of the ntuple size
  // Fill all generated particles
  bool do_mc_fill = false;
  if (!is_real_data && do_mc_fill){
    genps_p4 = tas::genps_p4();
    genps_id = tas::genps_id();
    genps_id_mother = tas::genps_id_mother();
    genps_idx_mother = tas::genps_idx_mother();
    genps_status = tas::genps_status();
    genps_id_grandma = tas::genps_id_simplegrandma();
    for (unsigned int i = 0; i < tas::els_p4().size(); i++){
        int mcidx = tas::els_mc3idx().at(i);
        if (mcidx < 0) continue;
        all_leps_p4.push_back(tas::els_p4()[i]);
        all_leps_id.push_back(-11*tas::els_charge().at(i));
        all_leps_mc3idx.push_back(mcidx);
        all_leps_veto.push_back(isGoodVetoElectron(i));
    }
    for (unsigned int i = 0; i < tas::mus_p4().size(); i++){
        int mcidx = tas::mus_mc3idx().at(i);
        if (mcidx < 0) continue;
        all_leps_p4.push_back(tas::mus_p4()[i]);
        all_leps_id.push_back(-13*tas::mus_charge().at(i));
        all_leps_mc3idx.push_back(mcidx);
        all_leps_veto.push_back(isGoodVetoMuon(i));
    }
  }

  if (!is_real_data) {
      // compute how many gen b-jets there are other than from top
      // (i.e., match status 1 gen b-quarks to genjets and count)
      extragenb = 0;
      ngenjets = 0;
      ngenjets30 = 0;
      genht = 0;
      genht30 = 0;
      std::vector<int> bidxs; // check for duplicates so we don't double count/match
      extragenbmoms.clear();
      for (unsigned int ij = 0; ij < tas::genjets_p4NoMuNoNu().size(); ij++){
          ngenjets++;
          auto jet = tas::genjets_p4NoMuNoNu()[ij];
          auto pt = jet.pt();
          bool foundB = false;
          float mindR = 0.25;
          int bestidx = -1;

          // Note, for wjets HT stitching, it's actually Events.Draw("Sum$(genps_p4.pt()*(genps_status==23)*(abs(genps_id)<=5 || genps_id==21))","","",10000)
          
          genht += pt;
          if (pt > 30.) {
              ngenjets30++;
              if (fabs(jet.eta()) < 2.4) {
                  genht30 += pt;
              }
          }

          for (unsigned int igen = 0; igen < tas::genps_p4().size(); igen++){
              auto genp4 = tas::genps_p4()[igen];
              int stat = tas::genps_status()[igen];
              int id = tas::genps_id()[igen];
              int mid = tas::genps_id_mother()[igen];
              if (stat == 1111) continue;
              // if (stat != 1) continue; // <-- literally nothing useful is status 1, so ignoring
              if (abs(id) != 5) continue;
              if (!idIsBeauty(abs(id))) continue;
              if (abs(mid) == 6) continue;
              if (abs(mid) == 22) continue;
              if (abs(mid) == 23) continue;
              if (abs(mid) == 24) continue;
              if (abs(mid) == 25) continue;

              if (std::find(bidxs.begin(), bidxs.end(), igen) != bidxs.end()) continue;
              float dR = ROOT::Math::VectorUtil::DeltaR(genp4,jet);
              if (dR < mindR) {
                  foundB = true;
                  bestidx = igen;
              }
          }
          if (foundB) {
              extragenb += 1;
              bidxs.push_back(bestidx);
              extragenbmoms.push_back(tas::genps_id_mother()[bestidx]);
          }
      }
  }

  if (!is_real_data) {
      // compute how many gen c-jets there are other than from top
      // (i.e., match status 1 gen c-quarks to genjets and count)
      extragenc = 0;
      std::vector<int> cidxs; // check for duplicates so we don't double count/match
      extragencmoms.clear();
      for (unsigned int ij = 0; ij < tas::genjets_p4NoMuNoNu().size(); ij++){
          auto jet = tas::genjets_p4NoMuNoNu()[ij];
          auto pt = jet.pt();
          bool foundC = false;
          float mindR = 0.25;
          int bestidx = -1;

          for (unsigned int igen = 0; igen < tas::genps_p4().size(); igen++){
              auto genp4 = tas::genps_p4()[igen];
              int stat = tas::genps_status()[igen];
              int id = tas::genps_id()[igen];
              int mid = tas::genps_id_mother()[igen];
              if (stat == 1111) continue;
              // if (stat != 1) continue; // <-- literally nothing useful is status 1, so ignoring
              if (abs(id) != 4) continue;
              if (!idIsCharm(abs(id))) continue;
              if (abs(mid) == 6) continue;
              if (abs(mid) == 22) continue;
              if (abs(mid) == 23) continue;
              if (abs(mid) == 24) continue;
              if (abs(mid) == 25) continue;

              if (std::find(cidxs.begin(), cidxs.end(), igen) != cidxs.end()) continue;
              float dR = ROOT::Math::VectorUtil::DeltaR(genp4,jet);
              if (dR < mindR) {
                  foundC = true;
                  bestidx = igen;
              }
          }
          if (foundC) {
              extragenc += 1;
              cidxs.push_back(bestidx);
              extragencmoms.push_back(tas::genps_id_mother()[bestidx]);
          }
      }
  }



  if (!is_real_data && iSignal > 0) {
      ncharginos = 0;
      nhiggs = 0;
      for (unsigned int gp=0;gp<tas::genps_id().size();gp++) {
          // for T5qqqqVV
          if (abs(tas::genps_id()[gp])==1000024 && tas::genps_status()[gp]==22) ncharginos++;
          // for T6ttHZ
          if (abs(tas::genps_id()[gp])==25 && tas::genps_status()[gp]==22) nhiggs++;
      }
  }

  //Determine and save jet and b-tag variables, raw
  std::pair <vector <Jet>, vector <Jet> > jet_results = SSJetsCalculator(jetCorr, 2, isFastsim);
  vector <LorentzVector> jets_raw;
  vector <LorentzVector> btags_raw;
  vector <float> jets_undoJEC_raw;
  for (unsigned int i = 0; i < jet_results.first.size(); i++)  jets_raw.push_back(jet_results.first.at(i).p4());
  for (unsigned int i = 0; i < jet_results.second.size(); i++) btags_raw.push_back(jet_results.second.at(i).p4());
  for (unsigned int i = 0; i < jet_results.first.size(); i++)  jets_undoJEC_raw.push_back(jet_results.first.at(i).undo_jec());
  njets_raw = jets_raw.size();
  nbtags_raw = btags_raw.size();
  ht_raw = 0;
  for (unsigned int i = 0; i < jets_raw.size(); i++) ht_raw += jets_raw.at(i).pt()*jets_undoJEC_raw.at(i);


  //Determine and save jet and b-tag variables, corrected
  jet_results = SSJetsCalculator(jetCorr, 1);
  for (unsigned int i = 0; i < jet_results.first.size(); i++) jets.push_back(jet_results.first.at(i).p4());
  for (unsigned int i = 0; i < jet_results.first.size(); i++) jets_flavor.push_back(jet_results.first.at(i).mcFlavor()); // NJA
  for (unsigned int i = 0; i < jet_results.first.size(); i++) jets_disc.push_back(jet_results.first.at(i).disc());
  vector<int> good_ijqs;
  vector<int> good_ijbs;
  for (unsigned int i = 0; i < jet_results.first.size(); i++) {
      auto jet = jet_results.first.at(i);
      if (!jet.isBtag()) good_ijqs.push_back(jet.idx());
  }
  for (unsigned int i = 0; i < jet_results.second.size(); i++) {
      auto jet = jet_results.second.at(i);
      good_ijbs.push_back(jet.idx());
  }


  for (unsigned int i = 0; i < jet_results.first.size(); i++) jets_JEC.push_back(jet_results.first.at(i).jec());
  for (unsigned int i = 0; i < jet_results.first.size(); i++) jets_undoJEC.push_back(jet_results.first.at(i).undo_jec());

  // for (unsigned int i = 0; i < jet_results.first.size(); i++) {
  //     auto jet = jet_results.first.at(i);
  //     float bsf = 1.0;
  //     if (gconf.year != 2016) {
  //         float jet_pt = jet.p4().pt()*jet.undo_jec()*jet.jec();
  //         float jet_eta = jet.p4().eta();
  //         int jet_mcFlavour = jet.mcFlavor();
  //         BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG;
  //         if (abs(jet_mcFlavour) == 5) flavor = BTagEntry::FLAV_B;
  //         else if (abs(jet_mcFlavour) == 4) flavor = BTagEntry::FLAV_C;
  //         bsf = btcr4->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
  //     }
  //     jets_bsf.push_back(bsf);
  // }

  for (unsigned int i = 0; i < jet_results.second.size(); i++) {
      auto jet = jet_results.second.at(i).p4();
      auto disc = jet_results.second.at(i).disc();
      // auto cdisc = jet_results.second.at(i).cdisc();
      btags.push_back(jet);
      btags_flavor.push_back(jet_results.second.at(i).mcFlavor());
      btags_disc.push_back(disc);
      // btags_cdisc.push_back(cdisc);
      btags_JEC.push_back(jet_results.second.at(i).jec());
      btags_undoJEC.push_back(jet_results.second.at(i).undo_jec());
  }

  for (unsigned int i = 0; i < jet_results.first.size(); i++){
    jecUnc->setJetEta(jets[i].eta());
    jecUnc->setJetPt(jets[i].pt()*jets_undoJEC[i]*jets_JEC[i]);
    jets_unc.push_back(jecUnc->getUncertainty(true));
  }
  for (unsigned int i = 0; i < jet_results.second.size(); i++){
    jecUnc->setJetEta(btags[i].eta());
    jecUnc->setJetPt(btags[i].pt()*btags_undoJEC[i]*btags_JEC[i]);
    btags_unc.push_back(jecUnc->getUncertainty(true));
  }
  njets = jets.size();
  nbtags = btags.size();
  ht = 0;
  for (unsigned int i = 0; i < jets.size(); i++){
     // jet_pt.push_back(jets.at(i).pt()*jets_undoJEC.at(i)*jets_JEC.at(i));
     ht += jets.at(i).pt()*jets_undoJEC.at(i)*jets_JEC.at(i);
  }

  if (verbose) cout << "ht: " << ht << endl;


  //now look at jets for get the btag scale factor (need to save down to 25 GeV corrected)
  // for applying btagging SFs, using Method 1a from the twiki below:
  //   https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods#1a_Event_reweighting_using_scale
  //   https://twiki.cern.ch/twiki/pub/CMS/BTagSFMethods/Method1aExampleCode_CSVM.cc.txt
  // Moved from additive to multiplicative uncertainties and split heavy/light by stealing from
  // https://github.com/cmstas/MT2Analysis/blob/master/babymaker/ScanChain.cc#L2235
  jet_results = SSJetsCalculator(jetCorr, 1, 0, 1);
  for (int reader_id=1; reader_id<=4; reader_id++) { // Loop over each id file
    float btagprob_data = 1.;
    float btagprob_err_heavy_UP = 1.;
    float btagprob_err_heavy_DN = 1.;
    float btagprob_err_light_UP = 1.;
    float btagprob_err_light_DN = 1.;
    float btagprob_mc = 1.;
    for (unsigned int i = 0; i < jet_results.first.size(); i++) {
       if (is_real_data) continue;
       //get btag eff weights
       float jet_pt = jet_results.first.at(i).p4().pt()*jet_results.first.at(i).undo_jec()*jet_results.first.at(i).jec();
       // Don't consider any jets below 25
       if (jet_pt<bjetMinPt) continue;
       // Don't consider nonbjets with 25<pT<40
       // if ((!jet_results.first.at(i).isBtag()) && (jet_pt<40.)) continue;
       float jet_eta = jet_results.first.at(i).p4().eta();
       int jet_mcFlavour = jet_results.first.at(i).mcFlavor();
       float eff = getBtagEffFromFile(jet_pt, jet_eta, jet_mcFlavour);
       BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG;
       if (abs(jet_mcFlavour) == 5) flavor = BTagEntry::FLAV_B;
       else if (abs(jet_mcFlavour) == 4) flavor = BTagEntry::FLAV_C;
       // float pt_cutoff = std::max(30.,std::min(669.,double(jet_pt)));
       float weight_cent(1.), weight_UP(1.), weight_DN(1.);
       if (gconf.year == 2016) {
           if (gconf.cmssw_ver == 80) {
               if (reader_id == 1) {
                   weight_cent    = btcr1->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
                   weight_UP    = btcr1->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
                   weight_DN    = btcr1->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
               } else if (reader_id == 2) {
                   weight_cent    = btcr2->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
                   weight_UP    = btcr2->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
                   weight_DN    = btcr2->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
               }
           } else {
               if (reader_id == 1) {
                   weight_cent    = btcr1->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
                   weight_UP    = btcr1->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
                   weight_DN    = btcr1->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
               }
           }
       } else if (gconf.year == 2017) {
           if (reader_id == 1) {
               weight_cent    = btcr1->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
               weight_UP    = btcr1->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
               weight_DN    = btcr1->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
           } else if (reader_id == 2) {
               weight_cent    = btcr2->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
               weight_UP    = btcr2->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
               weight_DN    = btcr2->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
           } else if (reader_id == 3) {
               weight_cent    = btcr3->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
               weight_UP    = btcr3->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
               weight_DN    = btcr3->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
           } else if (reader_id == 4) {
               weight_cent    = btcr4->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
               weight_UP    = btcr4->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
               weight_DN    = btcr4->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
           }
       } else if (gconf.year == 2018) {
           if (reader_id == 1) {
               weight_cent    = btcr1->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
               weight_UP    = btcr1->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
               weight_DN    = btcr1->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
           }
       }
       if (isFastsim > 0) {
           weight_cent    = btcr_fs->eval_auto_bounds("central", flavor, jet_eta, jet_pt);
           weight_UP    = btcr_fs->eval_auto_bounds("up", flavor, jet_eta, jet_pt);
           weight_DN    = btcr_fs->eval_auto_bounds("down", flavor, jet_eta, jet_pt);
       }
       if (jet_results.first.at(i).isBtag()) {
         btagprob_data *= weight_cent * eff;
         btagprob_mc *= eff;
         if (flavor == BTagEntry::FLAV_UDSG) {
             btagprob_err_light_UP *= weight_UP * eff;
             btagprob_err_light_DN *= weight_DN * eff;
             btagprob_err_heavy_UP *= weight_cent * eff;
             btagprob_err_heavy_DN *= weight_cent * eff;
         }
         else {
             btagprob_err_light_UP *= weight_cent * eff;
             btagprob_err_light_DN *= weight_cent * eff;
             btagprob_err_heavy_UP *= weight_UP * eff;
             btagprob_err_heavy_DN *= weight_DN * eff;
         }
       }
       else {
         btagprob_data *= (1. - weight_cent * eff);
         btagprob_mc *= (1. - eff);
         if (flavor == BTagEntry::FLAV_UDSG) {
           btagprob_err_light_UP *= (1. - weight_UP * eff);
           btagprob_err_light_DN *= (1. - weight_DN * eff);
           btagprob_err_heavy_UP *= (1. - weight_cent * eff);
           btagprob_err_heavy_DN *= (1. - weight_cent * eff);
         }
         else {
           btagprob_err_light_UP *= (1. - weight_cent * eff);
           btagprob_err_light_DN *= (1. - weight_cent * eff);
           btagprob_err_heavy_UP *= (1. - weight_UP * eff);
           btagprob_err_heavy_DN *= (1. - weight_DN * eff);
         }
       }

       if (verbose) {
           cout << Form("proc jet pt, eta, isBtagged, mcFlav = %f, %f, %i, %i",jet_pt,jet_eta,jet_results.first.at(i).isBtag(),flavor) << endl;
           cout << Form("eff, SF = %f %f",eff,weight_cent) << endl;
           cout << Form("partial SF is %f",btagprob_data / btagprob_mc) << endl;
       }
    }
    if (gconf.year == 2016) {
        if (gconf.cmssw_ver == 80) {
        if (reader_id == 1) {
            weight_btagsf1 = btagprob_data / btagprob_mc;
            weight_btagsf1_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
            weight_btagsf1_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
            weight_btagsf1_light_UP = btagprob_err_light_UP / btagprob_mc;
            weight_btagsf1_light_DN = btagprob_err_light_DN / btagprob_mc;
        } else if (reader_id == 2) {
            weight_btagsf2 = btagprob_data / btagprob_mc;
            weight_btagsf2_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
            weight_btagsf2_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
            weight_btagsf2_light_UP = btagprob_err_light_UP / btagprob_mc;
            weight_btagsf2_light_DN = btagprob_err_light_DN / btagprob_mc;
        }
        } else {
            if (reader_id == 1) {
                weight_btagsf1 = btagprob_data / btagprob_mc;
                weight_btagsf1_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
                weight_btagsf1_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
                weight_btagsf1_light_UP = btagprob_err_light_UP / btagprob_mc;
                weight_btagsf1_light_DN = btagprob_err_light_DN / btagprob_mc;
            }
        }
    } else if (gconf.year == 2017) {
        if (reader_id == 1) {
            weight_btagsf1 = btagprob_data / btagprob_mc;
            weight_btagsf1_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
            weight_btagsf1_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
            weight_btagsf1_light_UP = btagprob_err_light_UP / btagprob_mc;
            weight_btagsf1_light_DN = btagprob_err_light_DN / btagprob_mc;
        } else if (reader_id == 2) {
            weight_btagsf2 = btagprob_data / btagprob_mc;
            weight_btagsf2_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
            weight_btagsf2_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
            weight_btagsf2_light_UP = btagprob_err_light_UP / btagprob_mc;
            weight_btagsf2_light_DN = btagprob_err_light_DN / btagprob_mc;
        } else if (reader_id == 3) {
            weight_btagsf3 = btagprob_data / btagprob_mc;
            weight_btagsf3_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
            weight_btagsf3_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
            weight_btagsf3_light_UP = btagprob_err_light_UP / btagprob_mc;
            weight_btagsf3_light_DN = btagprob_err_light_DN / btagprob_mc;
        } else if (reader_id == 4) {
            weight_btagsf4 = btagprob_data / btagprob_mc;
            weight_btagsf4_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
            weight_btagsf4_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
            weight_btagsf4_light_UP = btagprob_err_light_UP / btagprob_mc;
            weight_btagsf4_light_DN = btagprob_err_light_DN / btagprob_mc;
        }
    } else if (gconf.year == 2018) {
        if (reader_id == 1) {
            weight_btagsf1 = btagprob_data / btagprob_mc;
            weight_btagsf1_heavy_UP = btagprob_err_heavy_UP / btagprob_mc;
            weight_btagsf1_heavy_DN = btagprob_err_heavy_DN / btagprob_mc;
            weight_btagsf1_light_UP = btagprob_err_light_UP / btagprob_mc;
            weight_btagsf1_light_DN = btagprob_err_light_DN / btagprob_mc;
        }
    }
  }

  // Now we do this stupid thing of rolling a random number and
  // deciding which era the MC "belongs to", to assign the proper btag weight

  TRandom *tr1 = new TRandom(tas::evt_event());
  float rand = tr1->Rndm();

  if (gconf.year == 2016) {
      if (gconf.cmssw_ver == 80) {
          if (rand < 0.55) {
              // B-F is 55% of the lumi
              weight_btagsf =    weight_btagsf1;
              weight_btagsf_heavy_UP = weight_btagsf1_heavy_UP;
              weight_btagsf_heavy_DN = weight_btagsf1_heavy_DN;
              weight_btagsf_light_UP = weight_btagsf1_light_UP;
              weight_btagsf_light_DN = weight_btagsf1_light_DN;
          } else {
              weight_btagsf =    weight_btagsf2;
              weight_btagsf_heavy_UP = weight_btagsf2_heavy_UP;
              weight_btagsf_heavy_DN = weight_btagsf2_heavy_DN;
              weight_btagsf_light_UP = weight_btagsf2_light_UP;
              weight_btagsf_light_DN = weight_btagsf2_light_DN;
          }
      } else {
          weight_btagsf =    weight_btagsf1;
          weight_btagsf_heavy_UP = weight_btagsf1_heavy_UP;
          weight_btagsf_heavy_DN = weight_btagsf1_heavy_DN;
          weight_btagsf_light_UP = weight_btagsf1_light_UP;
          weight_btagsf_light_DN = weight_btagsf1_light_DN;
      }
  }

  if (gconf.year == 2017) {
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2017Analysis
      // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation94X
      // RunB 4.767
      // RunC 9.583
      // RunD 4.224
      // RunE 9.261 | 7.793fb up to run 304671 | 1.468 after
      // RunF 13.463
      // 3 regions: B: 4.767, C-halfE: 21.6, halfE-F: 14.931 --> fractions: < 0.1154 < 0.6385 < 1.0000
      // Seed by event number to keep deterministic
      // weight_btagsf (no number) is the chosen one based on the era
      // 1,2,3,4 are period1,period2,period3, or inclusive
      if (rand < 0.1154) {
          weight_btagsf =    weight_btagsf1;
          weight_btagsf_heavy_UP = weight_btagsf1_heavy_UP;
          weight_btagsf_heavy_DN = weight_btagsf1_heavy_DN;
          weight_btagsf_light_UP = weight_btagsf1_light_UP;
          weight_btagsf_light_DN = weight_btagsf1_light_DN;
      } else if (rand < 0.6385) {
          weight_btagsf =    weight_btagsf2;
          weight_btagsf_heavy_UP = weight_btagsf2_heavy_UP;
          weight_btagsf_heavy_DN = weight_btagsf2_heavy_DN;
          weight_btagsf_light_UP = weight_btagsf2_light_UP;
          weight_btagsf_light_DN = weight_btagsf2_light_DN;
      } else {
          weight_btagsf =    weight_btagsf3;
          weight_btagsf_heavy_UP = weight_btagsf3_heavy_UP;
          weight_btagsf_heavy_DN = weight_btagsf3_heavy_DN;
          weight_btagsf_light_UP = weight_btagsf3_light_UP;
          weight_btagsf_light_DN = weight_btagsf3_light_DN;
      }
  }

  if (gconf.year == 2018) {
      weight_btagsf =    weight_btagsf1;
      weight_btagsf_heavy_UP = weight_btagsf1_heavy_UP;
      weight_btagsf_heavy_DN = weight_btagsf1_heavy_DN;
      weight_btagsf_light_UP = weight_btagsf1_light_UP;
      weight_btagsf_light_DN = weight_btagsf1_light_DN;
  }

  if (false) {
      weight_btagsf_iter_central = 1.;
      weight_btagsf_iter_down_cferr1 = 1.;
      weight_btagsf_iter_down_cferr2 = 1.;
      weight_btagsf_iter_down_hf = 1.;
      weight_btagsf_iter_down_hfstats1 = 1.;
      weight_btagsf_iter_down_hfstats2 = 1.;
      weight_btagsf_iter_down_jes = 1.;
      weight_btagsf_iter_down_lf = 1.;
      weight_btagsf_iter_down_lfstats1 = 1.;
      weight_btagsf_iter_down_lfstats2 = 1.;
      weight_btagsf_iter_up_cferr1 = 1.;
      weight_btagsf_iter_up_cferr2 = 1.;
      weight_btagsf_iter_up_hf = 1.;
      weight_btagsf_iter_up_hfstats1 = 1.;
      weight_btagsf_iter_up_hfstats2 = 1.;
      weight_btagsf_iter_up_jes = 1.;
      weight_btagsf_iter_up_lf = 1.;
      weight_btagsf_iter_up_lfstats1 = 1.;
      weight_btagsf_iter_up_lfstats2 = 1.;
      for (unsigned int i = 0; i < jet_results.first.size(); i++) {
          auto jet = jet_results.first.at(i);
          if (is_real_data) continue;
          float jet_pt = jet.p4().pt()*jet.undo_jec()*jet.jec();
          if (jet_pt<bjetMinPt) continue;
          float jet_eta = jet.eta();
          float disc = jet.disc();
          int jet_mcFlavour = cms3.pfjets_hadronFlavour().at(jet.idx());
          BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG; // 2
          if (abs(jet_mcFlavour) == 5) flavor = BTagEntry::FLAV_B; // 0
          else if (abs(jet_mcFlavour) == 4) flavor = BTagEntry::FLAV_C; // 1
          weight_btagsf_iter_central *= btcr_test->eval_auto_bounds("central", flavor, jet_eta, jet_pt, disc);
          if (flavor == BTagEntry::FLAV_B) {
              weight_btagsf_iter_down_lf *= btcr_test->eval_auto_bounds("down_lf", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_up_lf *= btcr_test->eval_auto_bounds("up_lf", flavor, jet_eta, jet_pt, disc);
              // weight_btagsf_iter_down_hfstats1 *= btcr_test->eval_auto_bounds("down_hfstats1", flavor, jet_eta, jet_pt, disc);
              // weight_btagsf_iter_down_hfstats2 *= btcr_test->eval_auto_bounds("down_hfstats2", flavor, jet_eta, jet_pt, disc);
              // weight_btagsf_iter_up_hfstats1 *= btcr_test->eval_auto_bounds("up_hfstats1", flavor, jet_eta, jet_pt, disc);
              // weight_btagsf_iter_up_hfstats2 *= btcr_test->eval_auto_bounds("up_hfstats2", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_down_jes *= btcr_test->eval_auto_bounds("down_jes", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_up_jes *= btcr_test->eval_auto_bounds("up_jes", flavor, jet_eta, jet_pt, disc);
          }
          else if (flavor == BTagEntry::FLAV_C) {
              weight_btagsf_iter_down_cferr1 *= btcr_test->eval_auto_bounds("down_cferr1", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_down_cferr2 *= btcr_test->eval_auto_bounds("down_cferr2", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_up_cferr1 *= btcr_test->eval_auto_bounds("up_cferr1", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_up_cferr2 *= btcr_test->eval_auto_bounds("up_cferr2", flavor, jet_eta, jet_pt, disc);
          }
          else if (flavor == BTagEntry::FLAV_UDSG) {
              // weight_btagsf_iter_down_lfstats1 *= btcr_test->eval_auto_bounds("down_lfstats1", flavor, jet_eta, jet_pt, disc);
              // weight_btagsf_iter_down_lfstats2 *= btcr_test->eval_auto_bounds("down_lfstats2", flavor, jet_eta, jet_pt, disc);
              // weight_btagsf_iter_up_lfstats1 *= btcr_test->eval_auto_bounds("up_lfstats1", flavor, jet_eta, jet_pt, disc);
              // weight_btagsf_iter_up_lfstats2 *= btcr_test->eval_auto_bounds("up_lfstats2", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_down_hf *= btcr_test->eval_auto_bounds("down_hf", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_up_hf *= btcr_test->eval_auto_bounds("up_hf", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_down_jes *= btcr_test->eval_auto_bounds("down_jes", flavor, jet_eta, jet_pt, disc);
              weight_btagsf_iter_up_jes *= btcr_test->eval_auto_bounds("up_jes", flavor, jet_eta, jet_pt, disc);
          }
      }
  }

  // ISR stuff for 2016
  // https://github.com/cmstas/MT2Analysis/blob/cmssw80x/babymaker/ScanChain.cc#L1870
  // https://github.com/cmstas/CORE/blob/master/MCSelections.cc#L438
  if (!is_real_data) {
      nisrMatch = get_nisrMatch(jets);
      weight_isr = get_isrWeight(nisrMatch);
      float unc_isr = get_isrUnc(nisrMatch);
      weight_isr_UP = weight_isr + unc_isr;
      weight_isr_DN = weight_isr - unc_isr;
      weight_isr_dy = 1.00; weight_isr_tt = 1.00;
      if (nisrMatch == 0) { weight_isr_dy = 1.00; weight_isr_tt = 1.00; }
      if (nisrMatch == 1) { weight_isr_dy = 0.96; weight_isr_tt = 0.86; }
      if (nisrMatch == 2) { weight_isr_dy = 0.96; weight_isr_tt = 0.79; }
      if (nisrMatch == 3) { weight_isr_dy = 0.95; weight_isr_tt = 0.75; }
      if (nisrMatch == 4) { weight_isr_dy = 1.05; weight_isr_tt = 0.77; }
      if (nisrMatch == 5) { weight_isr_dy = 1.11; weight_isr_tt = 0.74; }
      if (nisrMatch == 6) { weight_isr_dy = 1.22; weight_isr_tt = 0.73; }
      if (nisrMatch == 7) { weight_isr_dy = 1.34; weight_isr_tt = 0.56; }
  }

  //Save Most jets
  vector <Jet> mostJets_jet;
  vector <int> mostJets_idx;
  vector <Jet> mostJets_jet_up;
  vector <Jet> mostJets_jet_dn;
  for (unsigned int i = 0; i < tas::pfjets_p4().size(); i++){
    //Alias
    LorentzVector jet = tas::pfjets_p4().at(i);

    //Cuts
    // if (jet.pt() < 5.) continue;
    if (fabs(jet.eta()) > 2.4) continue;

    //Calculate raw pt
    float rawpt = jet.pt()*tas::pfjets_undoJEC().at(i);

    //Calculate jet corr
    jetCorr->setJetEta(jet.eta());
    jetCorr->setJetPt(rawpt);
    jetCorr->setJetA(tas::pfjets_area().at(i));
    jetCorr->setRho(tas::evt_fixgridfastjet_all_rho());
    float JEC = jetCorr->getCorrection();

    //Calculate jet unc
    jecUnc->setJetEta(pfjets_p4().at(i).eta());
    jecUnc->setJetPt(rawpt*JEC);
    float jetUnc = jecUnc->getUncertainty(true);

    if ( (rawpt*JEC*(1-jetUnc) < bjetMinPt)
      && (rawpt*JEC*(1+jetUnc) < bjetMinPt)
      && (rawpt*JEC*(1       ) < bjetMinPt) ) continue;

    //Save results
    mostJets.push_back(jet);
    // if (rawpt*JEC > jetMinPt && pass_SS_jetID(i,isFastsim))            mostJets_jet.push_back(Jet(i, JEC));
    if (rawpt*JEC > bjetMinPt && pass_SS_jetID(i,isFastsim))            mostJets_jet.push_back(Jet(i, JEC));
    else                                                     mostJets_jet.push_back(Jet(-1, -9999));
    if (rawpt*JEC*(1+jetUnc) > bjetMinPt && pass_SS_jetID(i,isFastsim)) mostJets_jet_up.push_back(Jet(i, JEC));
    else                                                     mostJets_jet_up.push_back(Jet(-1, -9999));
    if (rawpt*JEC*(1-jetUnc) > bjetMinPt && pass_SS_jetID(i,isFastsim)) mostJets_jet_dn.push_back(Jet(i, JEC));
    else                                                     mostJets_jet_dn.push_back(Jet(-1, -9999));
    mostJets_rawp4.push_back(jet*tas::pfjets_undoJEC().at(i));
    mostJets_idx.push_back(i);
    mostJets_disc.push_back(Jet(i,JEC).disc());
    mostJets_unc.push_back(jetUnc);
    mostJets_JEC.push_back(JEC);
    mostJets_undoJEC.push_back(tas::pfjets_undoJEC().at(i));
    mostJets_passID.push_back(pass_SS_jetID(i,isFastsim));
  }

  //Separate cleaning for all the permutations
  mostJets_passCleaning = cleanJets(mostJets_jet);
  vector <bool> mostJets_up_passCleaning = cleanJets(mostJets_jet_up);
  vector <bool> mostJets_dn_passCleaning = cleanJets(mostJets_jet_dn);

  //Unc up/down results
  njetsAG = 0; // is the same as njets
  nbtagsAG = 0; // is the same as nbtags
  float htb = 0;
  float htb_jec_up = 0;
  float htb_jec_dn = 0;
  float htb_jer_up = 0;
  float htb_jer_dn = 0;
  int nlb40 = 0;
  int nlb40_jec_up = 0;
  int nlb40_jec_dn = 0;
  int nlb40_jer_up = 0;
  int nlb40_jer_dn = 0;
  int ntb40 = 0;
  int ntb40_jec_up = 0;
  int ntb40_jec_dn = 0;
  int ntb40_jer_up = 0;
  int ntb40_jer_dn = 0;
  // XXX NOMINAL
  for (unsigned int i = 0; i < mostJets.size(); i++){
    if (mostJets_passCleaning.at(i) == 0) continue;
    if (mostJets_passID.at(i) == 0) continue;
    float jet_pt = mostJets.at(i).pt()*mostJets_undoJEC.at(i)*mostJets_JEC.at(i);
    if (jet_pt > bjetMinPt && mostJets_disc.at(i) > gconf.btag_disc_wp) nbtagsAG++;
    if (jet_pt > bjetMinPt && mostJets_disc.at(i) > gconf.btag_disc_wp) htb += jet_pt;
    if (jet_pt > 40. && mostJets_disc.at(i) > gconf.WP_DEEPCSV_LOOSE) nlb40++;
    if (jet_pt > 40. && mostJets_disc.at(i) > gconf.WP_DEEPCSV_TIGHT) ntb40++;
    if (jet_pt > jetMinPt) njetsAG++;
  }
  // XXX JEC UP
  for (unsigned int i = 0; i < mostJets.size(); i++){
    if (mostJets_up_passCleaning.at(i) == 0) continue;
    if (mostJets_passID.at(i) == 0) continue;
    float jet_pt_up = mostJets.at(i).pt()*mostJets_undoJEC.at(i)*mostJets_JEC.at(i)*(1.0+mostJets_unc.at(i));
    if (jet_pt_up > jetMinPt) {
        njets_unc_up++;
        // jets_jec_up.push_back(mostJets.at(i));
        // jets_disc_up.push_back(mostJets_disc.at(i));
        ht_unc_up += jet_pt_up;
    }
    if ((jet_pt_up > 40.) && (mostJets_disc.at(i) > gconf.WP_DEEPCSV_LOOSE)) nlb40_jec_up++;
    if ((jet_pt_up > 40.) && (mostJets_disc.at(i) > gconf.WP_DEEPCSV_TIGHT)) ntb40_jec_up++;
    if (mostJets_disc.at(i) < gconf.btag_disc_wp) continue;
    if (jet_pt_up > bjetMinPt) {
        nbtags_unc_up++;
        htb_jec_up += jet_pt_up;
        // bjets_jec_up.push_back(mostJets.at(i));
    }
  }
  // XXX JEC DOWN
  for (unsigned int i = 0; i < mostJets.size(); i++){
    if (mostJets_dn_passCleaning.at(i) == 0) continue;
    if (mostJets_passID.at(i) == 0) continue;
    float jet_pt_dn = mostJets.at(i).pt()*mostJets_undoJEC.at(i)*mostJets_JEC.at(i)*(1.0-mostJets_unc.at(i));
    if (jet_pt_dn > jetMinPt) {
        njets_unc_dn++;
        // jets_jec_dn.push_back(mostJets.at(i));
        // jets_disc_dn.push_back(mostJets_disc.at(i));
        ht_unc_dn += jet_pt_dn;
    }
    if ((jet_pt_dn > 40.) && (mostJets_disc.at(i) > gconf.WP_DEEPCSV_LOOSE)) nlb40_jec_dn++;
    if ((jet_pt_dn > 40.) && (mostJets_disc.at(i) > gconf.WP_DEEPCSV_TIGHT)) ntb40_jec_dn++;
    if (mostJets_disc.at(i) < gconf.btag_disc_wp) continue;
    if (jet_pt_dn > bjetMinPt) {
        nbtags_unc_dn++;
        htb_jec_dn += jet_pt_dn;
        // bjets_jec_dn.push_back(mostJets.at(i));
    }
  }


    auto tmpup = getT1CHSMET_fromMINIAOD(jetCorr, jecUnc, 1, false, use_cleaned_met);
    auto tmpdn = getT1CHSMET_fromMINIAOD(jetCorr, jecUnc, 0, false, use_cleaned_met);
    met_unc_up = tmpup.first;
    metPhi_unc_up = tmpup.second;
    met_unc_dn = tmpdn.first;
    metPhi_unc_dn = tmpdn.second;

    // JetResolution
    
    //Save Most jets - JER
    if (!is_real_data) {
        vector <LorentzVector> jer_mostJets;
        vector <Jet> jer_mostJets_jet;
        vector <int> jer_mostJets_idx;
        vector <float> jer_mostJets_disc;
        vector <float> jer_mostJets_JEC;
        vector <float> jer_mostJets_JER;
        vector <float> jer_mostJets_undoJEC;
        vector <bool> jer_mostJets_passID;
        vector <Double_t> GenJetPt;
        res.resetSeed(tas::evt_event()); // reset seed for JER smearing
        res.loadVariable("Rho", cms3.evt_fixgridfastjet_all_rho());
        for (auto& genjet : cms3.genjets_p4NoMuNoNu()) {
            GenJetPt.push_back(genjet.pt());
        }
        for (unsigned int i = 0; i < tas::pfjets_p4().size(); i++){
            //Alias
            LorentzVector jet = tas::pfjets_p4().at(i);

            //Cuts
            // if (jet.pt() < 5.) continue;
            if (fabs(jet.eta()) > 2.4) continue;

            //Calculate raw pt
            float rawpt = jet.pt()*tas::pfjets_undoJEC().at(i);

            //Calculate jet corr
            jetCorr->setJetEta(jet.eta());
            jetCorr->setJetPt(rawpt);
            jetCorr->setJetA(tas::pfjets_area().at(i));
            jetCorr->setRho(tas::evt_fixgridfastjet_all_rho());
            float JEC = jetCorr->getCorrection();

            res.loadVariable("JetEta", pfjets_p4().at(i).eta());
            res.loadVariable("JetPt", rawpt*JEC); // corrected pT
            auto smearing = res.smear(cms3.pfjets_p4().at(i), cms3.genjets_p4NoMuNoNu(), GenJetPt, 0);
            float ptJER = smearing[0]; // these are the actual pTs, not scale factors, so don't multiply
            float JER = ptJER/(rawpt*JEC); // JER wrt *corrected* jet

            if ( (ptJER    < bjetMinPt) ) continue;

            //Save results
            jer_mostJets.push_back(jet);
            if (ptJER > bjetMinPt && pass_SS_jetID(i,isFastsim))            jer_mostJets_jet.push_back(Jet(i, JEC));
            else                                                     jer_mostJets_jet.push_back(Jet(-1, -9999));
            jer_mostJets_idx.push_back(i);
            jer_mostJets_disc.push_back(Jet(i,JEC).disc());
            jer_mostJets_JEC.push_back(JEC);
            jer_mostJets_JER.push_back(JER);
            jer_mostJets_undoJEC.push_back(tas::pfjets_undoJEC().at(i));
            jer_mostJets_passID.push_back(pass_SS_jetID(i,isFastsim));
        }

        // Cleaning
        vector <bool> jer_mostJets_passCleaning = cleanJets(jer_mostJets_jet);

        //Unc up/down results
        njets_JER_up  = 0;
        nbtags_JER_up = 0;
        ht_JER_up     = 0;
        for (unsigned int i = 0; i < jer_mostJets.size(); i++){
            if (jer_mostJets_passCleaning.at(i) == 0) continue;
            if (jer_mostJets_passID.at(i) == 0) continue;
            float jet_pt = jer_mostJets.at(i).pt()*jer_mostJets_undoJEC.at(i)*jer_mostJets_JEC.at(i)*jer_mostJets_JER.at(i);
            if (jet_pt > bjetMinPt && jer_mostJets_disc.at(i) > gconf.btag_disc_wp) {
                nbtags_JER_up++;
                htb_jer_up += jet_pt;
            }
            if (jet_pt > jetMinPt) ht_JER_up += jet_pt;
            if (jet_pt > jetMinPt) njets_JER_up++;
            if ((jet_pt > 40.) && (jer_mostJets_disc.at(i) > gconf.WP_DEEPCSV_LOOSE)) nlb40_jer_up++;
            if ((jet_pt > 40.) && (jer_mostJets_disc.at(i) > gconf.WP_DEEPCSV_TIGHT)) ntb40_jer_up++;
        }
        njets_JER_dn = max(2*njets-njets_JER_up,0);
        nbtags_JER_dn = max(2*nbtags-nbtags_JER_up,0);
        ht_JER_dn = max(2.*ht-ht_JER_up,0.);
        nlb40_jer_dn = max(2*nlb40-nlb40_jer_up,0);
        ntb40_jer_dn = max(2*ntb40-ntb40_jer_up,0);
        htb_jer_dn = max(2.*htb-htb_jer_up,0.);

        LorentzVector jetp4_jer_up(0,0,0,0);
        LorentzVector jetp4_jer_dn(0,0,0,0);
        res.resetSeed(tas::evt_event()); // reset seed for JER smearing
        for(unsigned int iJet = 0; iJet < cms3.pfjets_p4().size(); iJet++){
            res.loadVariable("JetEta", pfjets_p4().at(iJet).eta());
            res.loadVariable("JetPt", pfjets_p4().at(iJet).pt());
            auto smearing = res.smear(cms3.pfjets_p4().at(iJet), cms3.genjets_p4NoMuNoNu(), GenJetPt, 0);
            float ptJER = smearing[0]; // these are the actual pTs, not scale factors
            float JER = ptJER/pfjets_p4().at(iJet).pt();
            jetp4_jer_up += cms3.pfjets_p4().at(iJet) * (1.-1./JER);
            jetp4_jer_dn += cms3.pfjets_p4().at(iJet) * (1.-JER);
        }
        float metx = met*cos(metPhi);
        float mety = met*sin(metPhi);
        float metx_up = metx + jetp4_jer_up.px();
        float mety_up = mety + jetp4_jer_up.py();
        float metx_dn = metx + jetp4_jer_dn.px();
        float mety_dn = mety + jetp4_jer_dn.py();
        met_JER_up = sqrt(metx_up*metx_up+mety_up*mety_up);
        metPhi_JER_up = atan2(mety_up,metx_up);
        met_JER_dn = sqrt(metx_dn*metx_dn+mety_dn*mety_dn);
        metPhi_JER_dn = atan2(mety_dn,metx_dn);
    }

  //Verbose for jets
  if (verbose){
    cout << "njets: " << njets << endl;
    cout << "nbtags: " <<  nbtags << endl;
    cout << "Printing Jets: " << endl;
    for (unsigned int i = 0; i < jets.size(); i++) cout << i << " " << jets[i].pt() << " " << jets[i].eta() << endl;
  }

  //MT variables
  mt    = MT(lep1_coneCorrPt, lep1_p4.phi(), met, metPhi);
  mt_l2 = MT(lep2_coneCorrPt, lep2_p4.phi(), met, metPhi);
  mtmin = mt > mt_l2 ? mt_l2 : mt;

  mt2   = -1;
  if (lep3_idx >= 0 && (abs(lep3_id) == 11 || abs(lep3_id) == 13)) {
      if ((lep1_id>0) == (lep3_id>0)) {
          // all the same sign? just put it in highest MT2 bin
          mt2 =  300.;
      } else {
          if (lep1_coneCorrPt > lep2_coneCorrPt) mt2 = MT2(met, metPhi, lep1_p4*lep1_coneCorrPt/lep1_p4.pt(), lep3_p4*lep3_coneCorrPt/lep3_p4.pt());
          else mt2 = MT2(met, metPhi, lep2_p4*lep2_coneCorrPt/lep2_p4.pt(), lep3_p4*lep3_coneCorrPt/lep3_p4.pt());
      }
  }

  //Ht and MET
  if (verbose) cout << "ht: " << ht << endl;
  if (verbose) cout << "met: " << met << endl;

  //Make sure one of our triggers fired
  bool isLL = lep1_p4.pt() < 25. and lep2_p4.pt() < 25;

  // https://indico.cern.ch/event/718554/contributions/3027981/attachments/1667626/2674497/leptontriggerreview.pdf#page=12
  if (gconf.year == 2016) {
      if (isFastsim == 0) {
          if (isRunH) {
              if (passHLTTrigger(triggerName("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"))  ||
                      passHLTTrigger(triggerName("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v")))   (triggers |= 1<<2);
              if (passHLTTrigger(triggerName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")))              (triggers |= 1<<3);
              if (passHLTTrigger(triggerName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")))            (triggers |= 1<<4);
          } else {
              if (passHLTTrigger(triggerName("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"))  ||
                      passHLTTrigger(triggerName("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v")))   (triggers |= 1<<2);
              if (passHLTTrigger(triggerName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")))              (triggers |= 1<<3);
              if (passHLTTrigger(triggerName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v")))            (triggers |= 1<<4);
          }
          if (passHLTTrigger(triggerName("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")))        (triggers |= 1<<6);
          if (passHLTTrigger(triggerName("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300_v")))      (triggers |= 1<<0);
          if (passHLTTrigger(triggerName("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v")))    (triggers |= 1<<5);
          if (passHLTTrigger(triggerName("HLT_DoubleMu8_Mass8_PFHT300_v")))                      (triggers |= 1<<7);
          // JetHT triggers
          if (passHLTTrigger(triggerName("HLT_AK8PFJet450")))                      (triggers |= 1<<8);
          if (passHLTTrigger(triggerName("HLT_PFJet450")))                         (triggers |= 1<<9);
      }
      fired_trigger = false;
      if (triggers != 0) {
          if (ht<300) { // First for FT
              if ( hyp_type==0 && ((triggers & 1<<3)==(1<<3) || (triggers & 1<<4)==(1<<4)) ) fired_trigger = true;
              if ( hyp_type==3 && (triggers & 1<<6)==(1<<6) ) fired_trigger = true;
              if ( (hyp_type==1 || hyp_type==2) && ((triggers & 1<<1)==(1<<1) || (triggers & 1<<2)==(1<<2)) ) fired_trigger = true;
          }
          else {
              if ( hyp_type==0 && (triggers & 1<<7)==(1<<7) ) fired_trigger = true;
              if ( hyp_type==3 && (triggers & 1<<5)==(1<<5) ) fired_trigger = true;
              if ( (hyp_type==1 || hyp_type==2) && (triggers & 1<<0)==(1<<0) ) fired_trigger = true;
          }
          if (isJetHT) {
              if ( (triggers & 1<<8)==(1<<8) ) fired_trigger = true;
              if ( (triggers & 1<<9)==(1<<9) ) fired_trigger = true;
              // for JetHT, just throw away the event if we didn't trigger...
              if (!fired_trigger) return babyErrorStruct;
          }
      }
      if (triggers != 0) { // Again for SS
          // Dilepton-only (iso)  triggers for HH and HL. For LL, dilepton+HT (with HT>400 threshold since HLT threshold moved up)
          if ( hyp_type==0 && ((triggers & 1<<3)==(1<<3) || (triggers & 1<<4)==(1<<4)) ) fired_trigger_ss = true;
          if ( hyp_type==3 && (triggers & 1<<6)==(1<<6) ) fired_trigger_ss = true;
          if ( (hyp_type==1 || hyp_type==2) && ((triggers & 1<<1)==(1<<1) || (triggers & 1<<2)==(1<<2)) ) fired_trigger_ss = true;
          if (isLL) {
              if ( hyp_type==0 && (triggers & 1<<7)==(1<<7) ) fired_trigger_ss = true;
              if ( hyp_type==3 && (triggers & 1<<5)==(1<<5) ) fired_trigger_ss = true;
              if ( (hyp_type==1 || hyp_type==2) && (triggers & 1<<0)==(1<<0) ) fired_trigger_ss = true;
          }
      }
      if (isFastsim > 0) {
          // There's no trigger emulation in 80X fastsim, but we'll just set the branch to 1 for convenience
          fired_trigger = true;
      }
  } else if (gconf.year == 2017) {
      if (isFastsim == 0) {
          if (passHLTTrigger(triggerName("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))  ||
                  passHLTTrigger(triggerName("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")))   (triggers |= 1<<2);
          if (passHLTTrigger(triggerName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v")))              (triggers |= 1<<3);
          if (passHLTTrigger(triggerName("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")))        (triggers |= 1<<6);

          if (passHLTTrigger(triggerName("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v")) ||
                  passHLTTrigger(triggerName("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v")))      (triggers |= 1<<0);
          if (passHLTTrigger(triggerName("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v")))    (triggers |= 1<<5);
          if (passHLTTrigger(triggerName("HLT_DoubleMu4_Mass8_DZ_PFHT350_v")))                      (triggers |= 1<<7);
      }
      fired_trigger = false;
      if (triggers != 0) {
          // Pure dilepton triggers for all HT
          if ( hyp_type==0 && ((triggers & 1<<3)==(1<<3) || (triggers & 1<<4)==(1<<4)) ) fired_trigger = true;
          if ( hyp_type==3 && (triggers & 1<<6)==(1<<6) ) fired_trigger = true;
          if ( (hyp_type==1 || hyp_type==2) && ((triggers & 1<<1)==(1<<1) || (triggers & 1<<2)==(1<<2)) ) fired_trigger = true;

          // Store noniso triggers in secondary variable for convenience (even though we can always do the annoying trigger bit math)
      }
      if (triggers != 0) { // Again for SS
          // For 2017, 2018 had been using isolated triggers, so copy the result; reset and recalculate if LL
          fired_trigger_ss = fired_trigger;
          if (isLL) {
              fired_trigger_ss = false;
              if ( hyp_type==0 && (triggers & 1<<7)==(1<<7) ) fired_trigger_ss = true;
              if ( hyp_type==3 && (triggers & 1<<5)==(1<<5) ) fired_trigger_ss = true;
              if ( (hyp_type==1 || hyp_type==2) && (triggers & 1<<0)==(1<<0) ) fired_trigger_ss = true;
          }
      }
      if (isFastsim > 0) {
          // There's no trigger emulation in 80X fastsim, but we'll just set the branch to 1 for convenience
          fired_trigger = true;
      }
  } else if (gconf.year == 2018) {
      if (isFastsim == 0) {
          if (passHLTTrigger(triggerName("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"))  ||
                  passHLTTrigger(triggerName("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")))   (triggers |= 1<<2);
          if (passHLTTrigger(triggerName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v")))              (triggers |= 1<<3);
          if (passHLTTrigger(triggerName("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")))        (triggers |= 1<<6);

          if (passHLTTrigger(triggerName("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v")))      (triggers |= 1<<0);
          if (passHLTTrigger(triggerName("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v")))    (triggers |= 1<<5);
          if (passHLTTrigger(triggerName("HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v")))                      (triggers |= 1<<7);
      }
      fired_trigger = false;
      if (triggers != 0) {
          // Pure dilepton triggers for all HT
          if ( hyp_type==0 && ((triggers & 1<<3)==(1<<3) || (triggers & 1<<4)==(1<<4)) ) fired_trigger = true;
          if ( hyp_type==3 && (triggers & 1<<6)==(1<<6) ) fired_trigger = true;
          if ( (hyp_type==1 || hyp_type==2) && ((triggers & 1<<1)==(1<<1) || (triggers & 1<<2)==(1<<2)) ) fired_trigger = true;
      }
      if (triggers != 0) { // Again for SS
          // For 2017, 2018 had been using isolated triggers, so copy the result; reset and recalculate if LL
          fired_trigger_ss = fired_trigger;
          if (isLL) {
              fired_trigger_ss = false;
              if ( hyp_type==0 && (triggers & 1<<7)==(1<<7) ) fired_trigger_ss = true;
              if ( hyp_type==3 && (triggers & 1<<5)==(1<<5) ) fired_trigger_ss = true;
              if ( (hyp_type==1 || hyp_type==2) && (triggers & 1<<0)==(1<<0) ) fired_trigger_ss = true;
          }
      }
  }

  //Number of good vertices
  for (unsigned int i = 0; i < tas::vtxs_ndof().size(); i++){
    if (!isGoodVertex(i)) continue;
    nGoodVertices++;
  }

  if (!is_real_data) {
      if (gconf.year == 2016) {
          std::tie(prefire2016_sf, prefire2016_sfup, prefire2016_sfdown) = getPrefireInfo(gconf.year); // or gconf.year eventually;
          prefire_sf = prefire2016_sf;
          prefire_sfup = prefire2016_sfup;
          prefire_sfdown = prefire2016_sfdown;
      }
      else if (gconf.year == 2017) {
          std::tie(prefire2017_sf, prefire2017_sfup, prefire2017_sfdown) = getPrefireInfo(gconf.year);
          prefire_sf = prefire2017_sf;
          prefire_sfup = prefire2017_sfup;
          prefire_sfdown = prefire2017_sfdown;
      }
      // std::tie(prefire2017ele_sf, prefire2017ele_sferr, prefire2017ele_njets) = getPrefireInfo(2017, false); // use ele/pho map, not jet
  }

  bool failsFastSimJetFilter = 0;
  if(isFastsim) {

      // std::tie(prefire2016_sf, prefire2016_sferr, prefire2016_njets) = getPrefireInfo(2016); // or gconf.year eventually;
      // std::tie(prefire2017_sf, prefire2017_sferr, prefire2017_njets) = getPrefireInfo(2017);
      // std::tie(prefire2017ele_sf, prefire2017ele_sferr, prefire2017ele_njets) = getPrefireInfo(2017, false); // use ele/pho map, not jet

      for (unsigned int i = 0; i < mostJets.size(); i++){
          if (mostJets_passCleaning.at(i) == 0) continue;
          if (mostJets_passID.at(i) == 0) continue;
          float jet_pt = mostJets.at(i).pt()*mostJets_undoJEC.at(i)*mostJets_JEC.at(i);
          int jetIdx = mostJets_idx.at(i);
          bool isBtag = mostJets_disc.at(i) > gconf.btag_disc_wp;
          if( isBtag && jet_pt < bjetMinPt) continue;
          if(!isBtag && jet_pt < jetMinPt) continue;
          if(isBadFastsimJet(jetIdx)) {
              failsFastSimJetFilter = 1;
              break;
          }
      }
  }

  // Doesn't include 
  // passes_met_filters = (isFastsim > 0) ? !failsFastSimJetFilter : passesMETfiltersRun2(is_real_data);
  passes_met_filters = passesMETfiltersRun2(is_real_data, (isFastsim > 0));
  if (isFastsim > 0 and failsFastSimJetFilter) {
      passes_met_filters = false;
  }
  passes_chargedcand_filter = true;
  if (gconf.cmssw_ver == 80) {
      if (!badChargedCandidateFilterV2()) passes_chargedcand_filter = false;
  } else {
      if (!filt_BadChargedCandidateFilter()) passes_chargedcand_filter = false;
  }

  // if (gconf.year == 2016) {
  //     if (gconf.cmssw_ver == 94) {
  //         passes_met_filters = (isFastsim > 0) ? !failsFastSimJetFilter : passesMETfilters2016(is_real_data, /*ignorechargedcandfilter=*/true, /*ignorebadmuonfilter=*/true);
  //     } else {
  //         passes_met_filters = (isFastsim > 0) ? !failsFastSimJetFilter : passesMETfilters2016(is_real_data, /*ignorechargedcandfilter=*/true);
  //     }
  // } else if (gconf.year == 2017) {
  //     passes_met_filters = (isFastsim > 0) ? !failsFastSimJetFilter : passesMETfilters2017(is_real_data);
  // } else if (gconf.year == 2018) {
  //     passes_met_filters = (isFastsim > 0) ? !failsFastSimJetFilter : passesMETfilters2017(is_real_data);
  // }
  

  // if (is_real_data && isReMiniAOD) {
  //     evt_egclean_pfmet = tas::evt_egclean_pfmet();
  //     evt_muegclean_pfmet = tas::evt_muegclean_pfmet();
  //     evt_muegcleanfix_pfmet = tas::evt_muegcleanfix_pfmet();
  //     evt_uncorr_pfmet = tas::evt_uncorr_pfmet();
  //     filt_noBadMuons = tas::filt_noBadMuons();
  //     filt_duplicateMuons = tas::filt_duplicateMuons();
  //     filt_badMuons = tas::filt_badMuons();
  // }

  if (!is_real_data and (filename.find("TTTT") != std::string::npos)) {
      auto totp4 = LorentzVector(0,0,0,0);
      for (unsigned int igen = 0; igen < tas::genps_id().size(); igen++){
          int id = tas::genps_id()[igen];
          if (abs(id) != 6) continue;
          int stat = tas::genps_status()[igen];
          if (stat != 22) continue;
          totp4 += tas::genps_p4()[igen];
      }
      mfourtop = totp4.mass();
  }

  if (!is_real_data) {

      float pdgleptW = 0.3258;
      float genleptW = 1.0/3;
      nleptonicW = 0;
      nW = 0;
      nhadronicW = 0;
      for (unsigned int igen = 0; igen < tas::genps_p4().size(); igen++){
          int id = tas::genps_id()[igen];
          if (abs(id) != 24) continue;
          int smid = tas::genps_id_simplemother()[igen];
          if (abs(smid) == 24) continue;
          int stat = tas::genps_status()[igen];
          if (stat != 22 && stat != 52) continue;
          nW++;
      }
      // Sometimes, we lose a W in the event record for ttW. Pin it to 3 if it's not 3.
      if (nW != 3 && (filename.find("TTWJetsToLNu") != std::string::npos)) nW = 3;
      for (unsigned int igen = 0; igen < tas::genps_p4().size(); igen++){
          int id = tas::genps_id()[igen];
          int mid = tas::genps_id_mother()[igen];
          if (abs(id) != 12 && abs(id) != 14 && abs(id) != 16) continue;
          if (abs(mid) != 24) continue;
          if (!(tas::genps_isPromptFinalState()[igen])) continue;
          nleptonicW++;
      }
      nhadronicW = nW - nleptonicW;
      leptonicWSF = pow((pdgleptW/genleptW),nleptonicW);
      hadronicWSF = pow(((1-pdgleptW)/(1-genleptW)),nhadronicW);
      decayWSF = leptonicWSF*hadronicWSF;

      // Find 4 bquarks
      std::vector<int> genbidxs;
      for (unsigned int igen = 0; igen < tas::genps_p4().size(); igen++){
          int stat = tas::genps_status()[igen];
          int id = tas::genps_id()[igen];
          int mid = tas::genps_id_mother()[igen];
          if (abs(id) != 5) continue;
          if (abs(mid) != 6) continue;
          if (stat != 23) continue;
          genbidxs.push_back(igen);
      }
      // Find 4 pairs of daughters from W decays
      std::vector<std::pair<int,int> > genwdaughteridxs;
      std::vector<bool> genwleptonic;
      std::vector<bool> genwtau;
      for (unsigned int igen = 0; igen < tas::genps_p4().size(); igen++){
          int id = tas::genps_id()[igen];
          int mid = tas::genps_id_mother()[igen];
          if (abs(id) != 24) continue;
          if (abs(mid) != 6) continue;
          std::pair<int,int> daughter_indices = {-1,-1};
          for (unsigned int idau = 0; idau < tas::genps_p4().size(); idau++){
              int midx = tas::genps_idx_mother()[idau];
              if (midx != igen) continue;
              // Fill first, then second, then break once both are in, since we can have duplicates :(
              if (daughter_indices.first == -1) daughter_indices.first = idau;
              else if (daughter_indices.second == -1) daughter_indices.second = idau;
              else break;
          }
          if (daughter_indices.first == -1 && daughter_indices.second == -1) continue;
          genwdaughteridxs.push_back(daughter_indices);
          genwleptonic.push_back(abs(tas::genps_id()[daughter_indices.first]) <= 16 && abs(tas::genps_id()[daughter_indices.first]) >= 11);
          genwtau.push_back(abs(tas::genps_id()[daughter_indices.first]) <= 16 && abs(tas::genps_id()[daughter_indices.first]) >= 15);
      }
      // Only good if we found all b's and all W daughters
      gengood = genbidxs.size()==4 && genwleptonic.size()==4;
      nleptonic = std::accumulate(genwleptonic.begin(),genwleptonic.end(),0);
      ntau = std::accumulate(genwtau.begin(),genwtau.end(),0);
      if (gengood) {

          // genjets_p4 = genjets_p4NoMuNoNu();

          // Make vector of [12-2*nleptonic] quarks to be matched
          std::vector<int> qidxtomatch;
          for (auto ib : genbidxs) qidxtomatch.push_back(ib);
          for (unsigned int gwd = 0; gwd < genwdaughteridxs.size(); gwd++) {
              if (genwleptonic[gwd]) continue;
              qidxtomatch.push_back(genwdaughteridxs[gwd].first);
              qidxtomatch.push_back(genwdaughteridxs[gwd].second);
          }

          ndrlt0p4 = 0;
          for (unsigned int iq = 0; iq < qidxtomatch.size(); iq++) {
              for (unsigned int jq = iq+1; jq < qidxtomatch.size(); jq++) {
                  float dR = ROOT::Math::VectorUtil::DeltaR(tas::genps_p4()[qidxtomatch[iq]],tas::genps_p4()[qidxtomatch[jq]]);
                  if (dR < 0.4)  {
                      ndrlt0p4++;
                  }
              }
          }

          // Match pfjets to qidxs, so below vector indexes qidxtomatch
          std::vector<int> alreadymatched;
          jet_type.clear();
          bool ignore_already_found = true;
          for (unsigned int ij = 0; ij < jets.size(); ij++){
              auto jet = jets[ij];
              float mindR = 0.40;
              int bestidx = -1;
              bool foundgenjet = false;
              for (unsigned int iq = 0; iq < qidxtomatch.size(); iq++) {
                  auto genp4 = tas::genps_p4()[qidxtomatch[iq]];
                  if (!ignore_already_found && std::find(alreadymatched.begin(), alreadymatched.end(), iq) != alreadymatched.end()) continue;
                  float dR = ROOT::Math::VectorUtil::DeltaR(genp4,jet);
                  if (dR < mindR) {
                      bestidx = iq;
                      mindR = dR;
                  }
              }
              alreadymatched.push_back(bestidx);
              for (unsigned int igj = 0; igj < genjets_p4NoMuNoNu().size(); igj++){
                  auto genjet = tas::genjets_p4NoMuNoNu()[igj];
                  float dR = ROOT::Math::VectorUtil::DeltaR(genjet,jet);
                  if (dR < 0.4) { foundgenjet = true; break; }
              }
              if (bestidx < 0 && !foundgenjet) jet_type.push_back(-2); // other
              else if (bestidx < 0 && foundgenjet) jet_type.push_back(-1); // other
              else if (bestidx < 4) jet_type.push_back(1); // first 4 are b's
              else jet_type.push_back(2); // q from W
          }
          alreadymatched.clear();
          bjet_type.clear();
          for (unsigned int ij = 0; ij < btags.size(); ij++){
              auto jet = btags[ij];
              float mindR = 0.45;
              int bestidx = -1;
              bool foundgenjet = false;
              for (unsigned int iq = 0; iq < qidxtomatch.size(); iq++) {
                  auto genp4 = tas::genps_p4()[qidxtomatch[iq]];
                  if (!ignore_already_found && std::find(alreadymatched.begin(), alreadymatched.end(), iq) != alreadymatched.end()) continue;
                  float dR = ROOT::Math::VectorUtil::DeltaR(genp4,jet);
                  if (dR < mindR) {
                      bestidx = iq;
                      mindR = dR;
                  }
              }
              alreadymatched.push_back(bestidx);
              for (unsigned int igj = 0; igj < genjets_p4NoMuNoNu().size(); igj++){
                  auto genjet = tas::genjets_p4NoMuNoNu()[igj];
                  float dR = ROOT::Math::VectorUtil::DeltaR(genjet,jet);
                  if (dR < 0.4) { foundgenjet = true; break; }
              }
              if (bestidx < 0 && !foundgenjet) bjet_type.push_back(-2); // other
              else if (bestidx < 0 && foundgenjet) bjet_type.push_back(-1); // other
              else if (bestidx < 4) bjet_type.push_back(1); // first 4 are b's
              else bjet_type.push_back(2); // q from W
          }

      }

  }

  passfilter = 1;
  if (!is_real_data and (gconf.year == 2016)) {
      // don't take inclusive sample if 
      // ngenjets30>=7 && nleptonicW==2 && genht30>500
      // and take ht500 sample instead
      if (filename.find("TTTo2L2Nu_HT500Njet7") != std::string::npos) {
          // MCM configuration: http://uaf-8.t2.ucsd.edu/~namin/dis/?query=%2FTTTo2L2Nu_HT500Njet7_TuneCUETP8M2T4_13TeV-powheg-pythia8%2FRunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1%2FMINIAODSIM%2Cthis&type=mcm
          passfilter = (ngenjets30 >= 7 && nleptonicW==2 && genht30>500);
          scale1fb *= 1.4960; //derived in samples.py -- check again after making babies
      }
      if (filename.find("TTToSemiLepton_HT500Njet9") != std::string::npos) {
          passfilter = (ngenjets30 >= 9 && nleptonicW==1 && genht30>500);
          scale1fb *= 2.3486; //derived in samples.py -- check again after making babies
      }
      if (filename.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") != std::string::npos) {
          passfilter = !(ngenjets30 >= 7 && nleptonicW==2 && genht30>500) && !(ngenjets30 >= 9 && nleptonicW==1 && genht30>500);
      }
  }

  // FIXME
  hasgammatoll = false;
  if (!is_real_data && (filename.find("TTJets_") != std::string::npos)) {
      // /hadoop/cms/store/group/snt/run2_mc2018//TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1_MINIAODSIM_CMS4_V10-02-04/merged_ntuple_31.root
      vector<LorentzVector> cands_p4;
      vector<int> cands_id;
      int mtype = -1;
      for (unsigned int igen = 0; igen < tas::genps_p4().size(); igen++) {
          auto genp4 = tas::genps_p4()[igen];
          int id = tas::genps_id()[igen];
          int mid = tas::genps_id_mother()[igen];
          int gmid = tas::genps_id_simplegrandma()[igen];
          bool lep_is_emu = (abs(id) == 11) or (abs(id) == 13);
          bool mom_is_emu = (abs(mid) == 11) or (abs(mid) == 13);
          bool mom_is_pho = (abs(mid) == 22);
          bool gma_is_emu = (abs(gmid) == 11) or (abs(gmid) == 13);
          if (not lep_is_emu) continue;
          if (mom_is_emu or (mom_is_pho and gma_is_emu)) {
              if (mom_is_emu) mtype = 1;
              cands_p4.push_back(genp4);
              cands_id.push_back(id);
          }
      }
      if (cands_id.size() == 2) {
          if (cands_id[0] == -cands_id[1]) {
              hasgammatoll = true;
              if (mtype == 1) gammatollmomemu = true;
              gammatoll1 = cands_p4[0];
              gammatoll2 = cands_p4[1];
              gammatolldr = ROOT::Math::VectorUtil::DeltaR(gammatoll1,gammatoll2);
              std::cout << std::endl;
              std::cout <<  " hasgammatoll: " << hasgammatoll <<  " gammatolldr: " << gammatolldr <<  std::endl;
              std::cout << std::endl;
          }
      }
  }





  // Compute some variables to make it easier for draw statements and skimming
  skim = (njets_unc_dn>=2 or njets_JER_dn>=2 or njets>=2 or njets_unc_up>=2 or njets_JER_up>=2) and \
         (met_unc_dn>=50 or met_JER_dn>=50 or met>=50 or met_unc_up>=50 or met_JER_up>=50) and \
         (hyp_class != 4 or is_real_data);
  skim_nomet = (njets_unc_dn>=2 or njets_JER_dn>=2 or njets>=2 or njets_unc_up>=2 or njets_JER_up>=2) and \
         (hyp_class != 4 or is_real_data);
  br = passes_baseline_ft(njets, nbtags, met, ht, lep1_id, lep2_id, lep1_coneCorrPt, lep2_coneCorrPt);
  nleps = (lep3_passes_id and lep3_coneCorrPt > 20) ? ((lep4_passes_id and lep4_coneCorrPt > 20) ? 4 : 3) : 2;
  sr = signal_region_ft(njets, nbtags, met, ht, -1, lep1_id, lep2_id, lep1_coneCorrPt, lep2_coneCorrPt, lep3_coneCorrPt, nleps, hyp_class==6);

  float mtnonz = mtmin;
  if (hyp_class==6) {
      float zmass23 = lep2_id == -lep3_id ? (lep2_p4+lep3_p4).mass() : -999.0;
      float zmass31 = lep3_id == -lep1_id ? (lep3_p4+lep1_p4).mass() : -999.0;
      if (fabs(zmass31 - 91.2) < fabs(zmass23 - 91.2)) {
          mtnonz = MT(lep2_coneCorrPt, lep2_p4.phi(), met, metPhi);
      } else {
          mtnonz = MT(lep1_coneCorrPt, lep1_p4.phi(), met, metPhi);
      }
  }
  ss_nleps = (lep3_passes_id) ? (lep4_passes_id ? 4 : 3) : 2;
  ss_region = analysis_category_ss(lep1_id, lep2_id, lep1_coneCorrPt, lep2_coneCorrPt, lep3_coneCorrPt, ss_nleps, ht, met);
  ss_br = passes_baseline_ss(njets, nbtags, met, ht, lep1_id, lep2_id, lep1_coneCorrPt, lep1_coneCorrPt);
  if (ss_region != LowMet and met < 50.) ss_br = 0;
  ss_sr = signal_region_ss(njets, nbtags, met, ht, mtmin, lep1_id, lep2_id, lep1_coneCorrPt, lep2_coneCorrPt, lep3_coneCorrPt, ss_nleps, hyp_class==6, mtnonz);

  // BDT stuff
  float mjoverpt = 0.;
  for (unsigned int ijet = 0; ijet < jets.size(); ijet++) {
      const auto& jet = jets[ijet]*(jets_undoJEC[ijet]*jets_JEC[ijet]);
      mjoverpt = std::max(mjoverpt, jet.M()/jet.pt());
  }
  bdt_nbtags = nbtags;
  bdt_njets = njets;
  bdt_met = met;
  bdt_ptl2 = lep2_coneCorrPt;
  bdt_nlb40 = nlb40;
  bdt_ntb40 = ntb40;
  bdt_nleps = nleps;
  bdt_htb = htb;
  bdt_ml1j1 = (njets >= 1) ? (jets[0]*jets_undoJEC[0]*jets_JEC[0] + lep1_p4).M() : 0;
  bdt_dphil1l2 = DeltaPhi(lep1_p4.phi(),lep2_p4.phi());
  bdt_maxmjoverpt = mjoverpt;
  bdt_detal1l2 = fabs(lep1_p4.eta()-lep2_p4.eta());
  bdt_q1 = (lep1_id > 0) ? -1 : 1;
  bdt_ptj1 = (njets > 0) ? jets[0].pt()*jets_undoJEC[0]*jets_JEC[0] : 0.;
  bdt_ptj6 = (njets > 5) ? jets[5].pt()*jets_undoJEC[5]*jets_JEC[5] : 0.;
  bdt_ptj7 = (njets > 6) ? jets[6].pt()*jets_undoJEC[6]*jets_JEC[6] : 0.;
  bdt_ptj8 = (njets > 7) ? jets[7].pt()*jets_undoJEC[7]*jets_JEC[7] : 0.;
  bdt_ptl1 = lep1_coneCorrPt;
  bdt_ptl3 = lep3_coneCorrPt;

  // For up and down variations

  bdt_jec_up_nbtags = nbtags_unc_up;
  bdt_jec_dn_nbtags = nbtags_unc_dn;
  bdt_jer_up_nbtags = nbtags_JER_up;
  bdt_jer_dn_nbtags = nbtags_JER_dn;

  bdt_jec_up_njets = njets_unc_up;
  bdt_jec_dn_njets = njets_unc_dn;
  bdt_jer_up_njets = njets_JER_up;
  bdt_jer_dn_njets = njets_JER_dn;

  bdt_jec_up_met = met_unc_up;
  bdt_jec_dn_met = met_unc_dn;
  bdt_jer_up_met = met_JER_up;
  bdt_jer_dn_met = met_JER_dn;

  bdt_jec_up_htb = htb_jec_up;
  bdt_jec_dn_htb = htb_jec_dn;
  bdt_jer_up_htb = htb_jer_up;
  bdt_jer_dn_htb = htb_jer_dn;

  bdt_jec_up_nlb40 = nlb40_jec_up;
  bdt_jec_dn_nlb40 = nlb40_jec_dn;
  bdt_jer_up_nlb40 = nlb40_jer_up;
  bdt_jer_dn_nlb40 = nlb40_jer_dn;

  bdt_jec_up_ntb40 = ntb40_jec_up;
  bdt_jec_dn_ntb40 = ntb40_jec_dn;
  bdt_jer_up_ntb40 = ntb40_jer_up;
  bdt_jer_dn_ntb40 = ntb40_jer_dn;

  bdt_disc = get_prediction( bdt_nbtags, bdt_njets, bdt_met, bdt_ptl2, bdt_nlb40, bdt_ntb40, bdt_nleps, bdt_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_disc_jec_up = get_prediction( bdt_jec_up_nbtags, bdt_jec_up_njets, bdt_jec_up_met, bdt_ptl2, bdt_jec_up_nlb40, bdt_jec_up_ntb40, bdt_nleps, bdt_jec_up_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_disc_jer_up = get_prediction( bdt_jer_up_nbtags, bdt_jer_up_njets, bdt_jer_up_met, bdt_ptl2, bdt_jer_up_nlb40, bdt_jer_up_ntb40, bdt_nleps, bdt_jer_up_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_disc_jec_dn = get_prediction( bdt_jec_dn_nbtags, bdt_jec_dn_njets, bdt_jec_dn_met, bdt_ptl2, bdt_jec_dn_nlb40, bdt_jec_dn_ntb40, bdt_nleps, bdt_jec_dn_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_disc_jer_dn = get_prediction( bdt_jer_dn_nbtags, bdt_jer_dn_njets, bdt_jer_dn_met, bdt_ptl2, bdt_jer_dn_nlb40, bdt_jer_dn_ntb40, bdt_nleps, bdt_jer_dn_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);

  bdt_run2_disc = bdtrun2::get_prediction( bdt_nbtags, bdt_njets, bdt_met, bdt_ptl2, bdt_nlb40, bdt_ntb40, bdt_nleps, bdt_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_run2_disc_jec_up = bdtrun2::get_prediction( bdt_jec_up_nbtags, bdt_jec_up_njets, bdt_jec_up_met, bdt_ptl2, bdt_jec_up_nlb40, bdt_jec_up_ntb40, bdt_nleps, bdt_jec_up_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_run2_disc_jer_up = bdtrun2::get_prediction( bdt_jer_up_nbtags, bdt_jer_up_njets, bdt_jer_up_met, bdt_ptl2, bdt_jer_up_nlb40, bdt_jer_up_ntb40, bdt_nleps, bdt_jer_up_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_run2_disc_jec_dn = bdtrun2::get_prediction( bdt_jec_dn_nbtags, bdt_jec_dn_njets, bdt_jec_dn_met, bdt_ptl2, bdt_jec_dn_nlb40, bdt_jec_dn_ntb40, bdt_nleps, bdt_jec_dn_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);
  bdt_run2_disc_jer_dn = bdtrun2::get_prediction( bdt_jer_dn_nbtags, bdt_jer_dn_njets, bdt_jer_dn_met, bdt_ptl2, bdt_jer_dn_nlb40, bdt_jer_dn_ntb40, bdt_nleps, bdt_jer_dn_htb, bdt_q1, bdt_ptj1, bdt_ptj6, bdt_ptj7, bdt_ml1j1, bdt_dphil1l2, bdt_maxmjoverpt, bdt_ptl1, bdt_detal1l2, bdt_ptj8, bdt_ptl3);

  bdt_run2_SRDISC = getBDTBin_17(bdt_run2_disc, hyp_class==6);

  yearlumi = getLumi(year);

  lep1_pt = lep1_p4.pt();
  lep2_pt = lep2_p4.pt();
  lep3_pt = lep3_p4.pt();

  mass13 = ((lep3_pt>5.) ? (lep3_p4+lep1_p4).mass() : -1.);
  mass23 = ((lep3_pt>5.) ? (lep2_p4+lep3_p4).mass() : -1.);
  mass123 = ((lep3_pt>5.) ? (lep1_p4+lep2_p4+lep3_p4).mass() : -1.);

  lep1_eta = lep1_p4.eta();
  lep2_eta = lep2_p4.eta();
  lep3_eta = lep3_p4.eta();

  lep1_phi = lep1_p4.phi();
  lep2_phi = lep2_p4.phi();
  lep3_phi = lep3_p4.phi();

  if (!is_real_data) {
      weight_lepsf1 = leptonScaleFactor(year, lep1_id, lep1_coneCorrPt, lep1_eta, ht, FTANA);
      weight_lepsf2 = leptonScaleFactor(year, lep2_id, lep2_coneCorrPt, lep2_eta, ht, FTANA);
      weight_lepsf3 = ((lep3_passes_id && lep3_coneCorrPt > 20.) ? leptonScaleFactor(year, lep3_id, lep3_pt, lep3_eta, ht, FTANA) : 1);
      weight_lepsf = weight_lepsf1 * weight_lepsf2 * weight_lepsf3;
      weight_triggersf = 1.; // FIXME
      // weight_triggersf = triggerScaleFactor(year, lep1_id, lep2_id, lep1_pt, lep2_pt, lep1_eta, lep2_eta, ht);
      weight_pu = getTruePUw(year, trueNumInt[0], 0);
      weight_isrsf = 1.;
      if (filename.find("TTWJetsToLNu") != std::string::npos) weight_isrsf = isrWeight(year, nisrMatch, 1);
      if (filename.find("TTZToLL") != std::string::npos) weight_isrsf = isrWeight(year, nisrMatch, 2);
      if (filename.find("TTJets") != std::string::npos) weight_isrsf = isrWeight(year, nisrMatch, 10);
      if (filename.find("TT_TuneCUETP8M2T4_13TeV-powheg") != std::string::npos) weight_isrsf = isrWeight(year, nisrMatch, 10);
      weight = yearlumi * scale1fb * weight_lepsf * weight_triggersf * weight_pu * weight_btagsf * decayWSF * weight_isrsf;
  } else {
      weight = 1;
  }

  // std::cout <<  " run: " << run <<  " lumi: " << lumi <<  " event: " << event <<  " lep1_id: " << lep1_id <<  " lep2_id: " << lep2_id <<  " lep1_pt: " << lep1_pt <<  " lep2_pt: " << lep2_pt <<  std::endl;

  //Fill Baby
  BabyTree->Fill();

  return babyErrorStruct;

}

float babyMaker::getBtagEffFromFile(float pt, float eta, int mcFlavour){
  if(!h_btag_eff_b || !h_btag_eff_c || !h_btag_eff_udsg) {
    std::cout << "babyMaker::getBtagEffFromFile: ERROR: missing input hists" << std::endl;
    return 1.;
  }

  // only use pt bins up to 400 GeV for charm and udsg
  float pt_cutoff = std::max(20.,std::min(399.,double(pt)));
  TH2D* h(0);
  if (abs(mcFlavour) == 5) {
    h = h_btag_eff_b;
    // use pt bins up to 600 GeV for b
    pt_cutoff = std::max(20.,std::min(599.,double(pt)));
  }
  else if (abs(mcFlavour) == 4) {
    h = h_btag_eff_c;
  }
  else {
    h = h_btag_eff_udsg;
  }

  int binx = h->GetXaxis()->FindBin(pt_cutoff);
  int biny = h->GetYaxis()->FindBin(fabs(eta));
  return h->GetBinContent(binx,biny);
}

