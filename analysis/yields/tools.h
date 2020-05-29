#ifndef ANALYSIS_YIELDS_TOOLS_H
#define ANALYSIS_YIELDS_TOOLS_H

//#include "TLorentzVector.h"
//#include "../misc/class_files/v8.02/SS.h"
#include "Math/LorentzVector.h"
#include <string>
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

/*
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
*/

float calcDeltaPhi(float phi1, float phi2);

float calcDeltaR(float eta1, float phi1, float eta2, float phi2);

float calcMT(float pt1, float phi1, float pt2, float phi2);

float minDR(const LorentzVector & lep , const vector<LorentzVector> & jet);

float PtMaxEta(const vector<LorentzVector> & jet);

struct Lepton
{
    Lepton () : name(""),v(LorentzVector(0,0,0,0)),id(0),isgood(false),isfo(false),ccpt(0.),miniiso(999.),dxy(999.),dz(999.),pt(0),eta(0),phi(0) {}
    Lepton (string name_, LorentzVector v_, int id_, bool isgood_, bool isfo_, float ccpt_, float miniiso_, float dxy_, float dz_) : 
            name(name_),v(v_),id(id_),isgood(isgood_),isfo(isfo_),ccpt(ccpt_),miniiso(miniiso_),dxy(dxy_),dz(dz_) {pt = v.pt(); eta = v.eta(); phi = v.phi();}
    string name;
    LorentzVector v;
    int id;
    bool isgood;
    bool isfo;
    float ccpt; 
    float miniiso;
    float dxy;
    float dz;   
    float pt;
    float eta;
    float phi;

    int print() {std::cout << this->name << ": pt = " << this->v.pt() << ", eta = " << this->v.eta() << ", is_good = " << this->isgood << ", isfo = " << this->isfo << std::endl;}
}; 

std::vector<std::pair<float, std::pair<int,int>>> MOSSF(vector<Lepton> lep);

std::pair<int,int> get_hyp_class(std::vector<Lepton> vleps, int ss_hyp_class);

bool is_ossf(Lepton lep1, Lepton lep2); 

bool sortByPt (float i, float j);

int printLeptons(std::vector<Lepton>);

/************************/
// signal regions
// 
// BRINCL = 1: SS or >=3, tight ID+SIO, pt > 25,20,20
//
//
// SSBR = 11: SS tight ID+ISO, pt > 25,20, veto third lepton passing FO selection
// SSBR_SF = 12: SS SF tight ID+ISO, pt > 25,20
// SSBR_OF = 13: OS SF tight ID+ISO, pt > 25,20
//
// SS_SR_NOB  = 21: SS TT Nb = 0 
// SS_SR_ONEB = 22: SS TT Nb = 1
// SS_SR_TWOB = 23: SS TT Nb >= 2
//
// SS_OS_NOB  = 31: OS TT Nb = 0 
// SS_OS_ONEB = 32: OS TT Nb = 1
// SS_OS_TWOB = 33: OS TT Nb >= 2
// 
// SS_TL_NOB  = 41: SS TL Nb = 0 
// SS_TL_ONEB = 42: SS TL Nb = 1
// SS_TL_TWOB = 43: SS TL Nb >= 2
//
// SS_LL_NOB  = 51: SS LL Nb = 0 
// SS_LL_ONEB = 52: SS LL Nb = 1
// SS_LL_TWOB = 53: SS LL Nb >= 2
//
//
// ML_BR  = 101: >=3L, tight ID+SIO, pt > 25,20,20
//
// ML_SR_NOZ_NOB   = 111: >=3L, Nz = 0, Nb = 0 
// ML_SR_NOZ_ONEB  = 112: >=3L, Nz = 0, Nb = 1 
// ML_SR_NOZ_TWOB   = 113: >=3L, Nz = 0, Nb = 2 
// ML_SR_ONEZ_NOB  = 114: >=3L, Nz = 1, Nb = 0 
// ML_SR_ONEZ_ONEB = 115: >=3L, Nz = 1, Nb = 1 
// ML_SR_ONEZ_TWOB = 116: >=3L, Nz = 1, Nb = 2
// ML_SR_TWOZ_NOB  = 117: >=3L, Nz = 2, Nb = 0 
// ML_SR_TWOZ_ONEB = 118: >=3L, Nz = 2, Nb = 1
// ML_SR_TWOZ_TWOB = 119: >=3L, Nz = 2, Nb = 2
//
// ML_TTL_NOZ_NOB   = 121: >=3L, Nz = 0, Nb = 0, TTL 
// ML_TTL_NOZ_ONEB  = 122: >=3L, Nz = 0, Nb = 1, TTL 
// ML_TTL_NOZ_TWOB   = 123: >=3L, Nz = 0, Nb = 2, TTL 
// ML_TTL_ONEZ_NOB  = 124: >=3L, Nz = 1, Nb = 0, TTL 
// ML_TTL_ONEZ_ONEB = 125: >=3L, Nz = 1, Nb = 1, TTL 
// ML_TTL_ONEZ_TWOB = 126: >=3L, Nz = 1, Nb = 2, TTL
// ML_TTL_TWOZ_NOB  = 127: >=3L, Nz = 2, Nb = 0, TTL 
// ML_TTL_TWOZ_ONEB = 128: >=3L, Nz = 2, Nb = 1, TTL
// ML_TTL_TWOZ_TWOB = 129: >=3L, Nz = 2, Nb = 2, TTL
//
// current (old) ML SRs
// ML_SR_NOB   = 131: >=3L, Nb = 0 
// ML_SR_ONEB  = 132: >=3L, Nb = 1 
// ML_SR_TWOB  = 133: >=3L, Nb = 2 
/************************/
/*
enum class Region {BRINCL = 1, SS_BR = 11, SS_BR_SF = 12, SS_BR_TL = 13, SS_SR_NOB = 21, SS_SR_ONEB = 22, SS_SR_TWOB = 23, SS_OS_NOB = 31, SS_OS_ONEB = 32, SS_OS_TWOB = 33, 
                    SS_LL_NOB = 41, SS_LL_ONEB = 42, SS_LL_TWOB = 43, SS_LL_NOB = 51, SS_LL_ONEB = 52, SS_LL_TWOB = 53, ML_BR = 101, 
                    ML_SR_NOZ_NOB   = 111,  ML_SR_NOZ_ONEB  = 112, ML_SR_NO_TWOB   = 113, ML_SR_ONEZ_NOB  = 114, ML_SR_ONEZ_ONEB = 115, ML_SR_ONEZ_TWOB = 116, ML_SR_TWOZ_NOB  = 117, ML_SR_TWOZ_ONEB = 118, ML_SR_TWOZ_TWOB = 119,
                    ML_TTL_NOZ_NOB   = 121,  ML_TTL_NOZ_ONEB  = 122, ML_TTL_NO_TWOB   = 123, ML_TTL_ONEZ_NOB  = 124, ML_TTL_ONEZ_ONEB = 125, ML_TTL_ONEZ_TWOB = 126, ML_TTL_TWOZ_NOB  = 127, ML_TTL_TWOZ_ONEB = 128, ML_TTL_TWOZ_TWOB = 129,
                    };

Region signal_region(int hyp_class, int nleps_tight, int nleps_fo, int nz, int nb, int nj) {

    z_index = std::min(nz,2);
    b_index = std::min(nb,2)+1;
    z_b_index = (nleps_fo > 2) ? 3*z_index + b_index : b_index;
    // ML first
    int ret = -1;
    if (nleps_tight >= 3) ret = 110;
    else if (nleps_tight == 2 and nleps_fo >= 3) ret = 120;
    else if (nleps_tight == 2 and hyp_class == 3) ret = 20;
    else if (hyp_class == 4) ret = 30;
    else if (hyp_class == 2) ret = 40;
    else if (hyp_class == 1) ret = 50;
    else return 1; 

    return ret+z_b_index;    
}
*/

// simpler function for region
int signal_region(int hyp_class, int nleps_tight, int nleps_fo, int nz, int nb, int nj); 

std::vector<std::pair<int,int>> num_z_cands(std::vector<LorentzVector> p4s, std::vector<int> ids);

#endif
