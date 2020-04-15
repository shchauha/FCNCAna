#ifndef _commonUtils_h_
#define _commonUtils_h_

// needed inside the below includes
enum ana_t { SSANA, FTANA };

#include "year_2016/all.h"
#include "year_2017/all.h"
#include "year_2018/all.h"
namespace yrun2 {
#include "year_run2/trigeffandsf.h"
}


float getTruePUw(int year, int nvtx, int which=0) { 
    // which = 0, 1, -1 for nominal, up, down
    if (which == 0) {
        if (year == 2016) return y2016::getTruePUw(nvtx);
        else if (year == 2017) return y2017::getTruePUw(nvtx);
        else if (year == 2018) return y2018::getTruePUw(nvtx);
    } else if (which == 1) {
        if (year == 2016) return y2016::getTruePUwUp(nvtx);
        else if (year == 2017) return y2017::getTruePUwUp(nvtx);
        else if (year == 2018) return y2018::getTruePUwUp(nvtx);
    } else if (which == -1) {
        if (year == 2016) return y2016::getTruePUwDn(nvtx);
        else if (year == 2017) return y2017::getTruePUwDn(nvtx);
        else if (year == 2018) return y2018::getTruePUwDn(nvtx);
    }
    return 0.;
}

float isrWeight(int year, int nisrmatch, int sample) { 
    if (year == 2016) return y2016::isrWeight(nisrmatch, sample);
    else if (year == 2017) return y2017::isrWeight(nisrmatch, sample);
    else if (year == 2018) return y2018::isrWeight(nisrmatch, sample);
    else return 0.;
}


float leptonScaleFactor(int year, int id, float pt, float eta, float ht) { 
    // FIXME Note: muon SFs do not go below pT of 20 - this is fine for four top, but
    // modify for SS!
    // Fix it and then remove protection lines in lepton_sf.h
    if (year == 2016) return y2016::leptonScaleFactor(id, pt, eta, ht);
    else if (year == 2017) return y2017::leptonScaleFactor(id, pt, eta, ht);
    else if (year == 2018) return y2018::leptonScaleFactor(id, pt, eta, ht);
    else return 0.;
}

float fastsim_leptonScaleFactor(int year, int id, float pt, float eta, float ht) { 
    // FIXME only 2016
    // modify for SS!
    if (year == 2016) return y2016::leptonScaleFactorFastSim(id, pt, eta, ht);
    else if (year == 2017) return y2017::leptonScaleFactorFastSim(id, pt, eta, ht);
    else if (year == 2018) return y2018::leptonScaleFactorFastSim(id, pt, eta, ht);
    // return y2016::leptonScaleFactorFastSim(id, pt, eta, ht);
    return 1.0;
}

float leptonScaleFactorError(int year, int id, float pt, float eta, float ht) { 
    if (year == 2016) return y2016::leptonScaleFactor_err(id, pt, eta, ht);
    else if (year == 2017) return y2017::leptonScaleFactor_err(id, pt, eta, ht);
    else if (year == 2018) return y2018::leptonScaleFactor_err(id, pt, eta, ht);
    else return 0.;
}

float triggerScaleFactor(int year, int pdgId1, int pdgId2, float pt1, float pt2, float eta1, float eta2, float ht, ana_t ana, int syst=0) { 
    // 0 - central, +1 and -1 for up and down variations
    // For 3 lepton events, don't scale, and take a 2% uncertainty
    if (ana == FTANA and year == 2016 and syst==0) return y2016::triggerScaleFactor(pdgId1, pdgId2, pt1, pt2, eta1, eta2, ht);
    if (syst == 2) syst = 0; // to bypass the above;
    return yrun2::TriggerWeight(pdgId1, pt1, eta1, pdgId2, pt2, eta2, ht, year, true, syst);
}

float fastsim_triggerScaleFactor(int year, int pdgId1, int pdgId2, float pt1, float pt2, float eta1, float eta2, float ht, int syst=0) {
    // FIXME only 2016 -- actually this is an efficiency since there is no trigger in 80X fastsim
    // modify for SS!
    // if (year == 2016) return y2016::FastSimTriggerEfficiency(id1,id2,pt1,pt2,eta1,eta2,ht);
    // return y2016::FastSimTriggerEfficiency(id1,id2,pt1,pt2,eta1,eta2,ht);
    return yrun2::TriggerWeight(pdgId1, pt1, eta1, pdgId2, pt2, eta2, ht, year, false, syst);
    // return 1.0;
}

float flipRate(int year, float pt, float eta) { 
    if (year == 2016) return y2016::flipRate(pt, eta);
    else if (year == 2017) return y2017::flipRate(pt, eta);
    else if (year == 2018) return y2018::flipRate(pt, eta);
    else return 0.;
}

float flipRateError(int year, float pt, float eta) { 
    if (year == 2016) return y2016::flipRateError(pt, eta);
    else if (year == 2017) return y2017::flipRateError(pt, eta);
    else if (year == 2018) return y2018::flipRateError(pt, eta);
    else return 0.;
}

float fakeRate(int year, int id, float pt, float eta, float ht, ana_t ana, bool newbins=false, bool allpt=false, bool isLL=false) { 
    if (year == 2016) return (newbins ? y2016::newbins::fakeRate(id, pt, eta, ht, ana, isLL) : y2016::fakeRate(id, pt, eta, ht, ana, isLL));
    else if (year == 2017) return (allpt ? y2017::allpt::fakeRate(id, pt, eta, ht, isLL) : y2017::fakeRate(id, pt, eta, ht, isLL));
    else if (year == 2018) return (allpt ? y2018::allpt::fakeRate(id, pt, eta, ht, isLL) : y2018::fakeRate(id, pt, eta, ht, isLL));
    else return 0.;
}
float fakeRateError(int year, int id, float pt, float eta, float ht, ana_t ana, bool newbins=false, bool allpt=false, bool isLL=false) { 
    if (year == 2016) return (newbins ? y2016::newbins::fakeRateError(id, pt, eta, ht, ana, isLL) : y2016::fakeRateError(id, pt, eta, ht, ana, isLL));
    else if (year == 2017) return (allpt ? y2017::allpt::fakeRateError(id, pt, eta, ht, isLL) : y2017::fakeRateError(id, pt, eta, ht, isLL));
    else if (year == 2018) return (allpt ? y2018::allpt::fakeRateError(id, pt, eta, ht, isLL) : y2018::fakeRateError(id, pt, eta, ht, isLL));
    else return 0.;
}

float alternativeFakeRate(int year, int id, float pt, float eta, float ht, ana_t ana, bool newbins=false, bool allpt=false, bool isLL=false) { 
    if (year == 2016) return (newbins ? y2016::newbins::alternativeFakeRate(id, pt, eta, ht, ana, isLL) : y2016::alternativeFakeRate(id, pt, eta, ht, ana, isLL));
    else if (year == 2017) return (allpt ? y2017::allpt::alternativeFakeRate(id, pt, eta, ht, isLL) : y2017::alternativeFakeRate(id, pt, eta, ht, isLL));
    else if (year == 2018) return (allpt ? y2018::allpt::alternativeFakeRate(id, pt, eta, ht, isLL) : y2018::alternativeFakeRate(id, pt, eta, ht, isLL));
    else return 0.;
}

float qcdMCFakeRate(int year, int id, float pt, float eta, float ht, ana_t ana, bool newbins=false, bool allpt=false, bool isLL=false) { 
    if (year == 2016) return (newbins ? y2016::newbins::qcdMCFakeRate(id, pt, eta, ht, ana, isLL) : y2016::qcdMCFakeRate(id, pt, eta, ht, ana, isLL));
    else if (year == 2017) return (allpt ? y2017::allpt::qcdMCFakeRate(id, pt, eta, ht, isLL) : y2017::qcdMCFakeRate(id, pt, eta, ht, isLL));
    else if (year == 2018) return (allpt ? y2018::allpt::qcdMCFakeRate(id, pt, eta, ht, isLL) : y2018::qcdMCFakeRate(id, pt, eta, ht, isLL));
    else return 0.;
}

float qcdMCFakeRateError(int year, int id, float pt, float eta, float ht, ana_t ana, bool newbins=false, bool allpt=false, bool isLL=false) { 
    if (year == 2016) return (newbins ? y2016::newbins::qcdMCFakeRateError(id, pt, eta, ht, ana, isLL) : y2016::qcdMCFakeRateError(id, pt, eta, ht, ana, isLL));
    else if (year == 2017) return (allpt ? y2017::allpt::qcdMCFakeRateError(id, pt, eta, ht, isLL) : y2017::qcdMCFakeRateError(id, pt, eta, ht, isLL));
    else if (year == 2018) return (allpt ? y2018::allpt::qcdMCFakeRateError(id, pt, eta, ht, isLL) : y2018::qcdMCFakeRateError(id, pt, eta, ht, isLL));
    else return 0.;
}

float getLumi(int year) {
    if (year == 2016) return y2016::getLumi();
    else if (year == 2017) return y2017::getLumi();
    else if (year == 2018) return y2018::getLumi();
    else return 0.;
}

#endif
