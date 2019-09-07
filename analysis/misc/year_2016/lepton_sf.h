#include "lepton_sfs_el.h"
#include "lepton_sfs_mu.h"

float leptonScaleFactor(int pdgId, float pt, float eta, float ht) {
  if (abs(pdgId)==13) return muonScaleFactor(pt, eta)*trackingSF(eta);
  else if (abs(pdgId)==11){
    if (ht>300) return electronScaleFactorHighHT(pt, eta)*electronGSF(pt, eta);
      else return electronScaleFactorLowHT(pt, eta)*electronGSF(pt, eta);
  }
  return 0.;
}

float leptonScaleFactor_err(int pdgId, float pt, float eta, float ht) {
  if (abs(pdgId)==13) return muonScaleFactor_error(pt, eta);
  else if (abs(pdgId)==11){
    if (ht>300) return electronScaleFactorHighHT_error(pt, eta);
      else return electronScaleFactorLowHT_error(pt, eta);
  }
  return 0.;
}

float leptonScaleFactor_legacy(int pdgId, float pt, float eta, float ht) {
  if (abs(pdgId)==13) return muonScaleFactor(pt, eta)*trackingSF(eta);
  else if (abs(pdgId)==11){
    return electronScaleFactor_legacy(pt, eta)*electronScaleFactorReco_legacy(pt, eta);
  }
  return 0.;
}

float leptonScaleFactor_err_legacy(int pdgId, float pt, float eta, float ht) {
    if (abs(pdgId)==13) return muonScaleFactor_error(pt, eta);
    else if (abs(pdgId)==11){
        float e1=electronScaleFactorError_legacy(pt,eta);
        float e2=electronScaleFactorRecoError_legacy(pt,eta);
        return pow(e1*e1+e2*e2,0.5);
    }
    return 0.;
}
