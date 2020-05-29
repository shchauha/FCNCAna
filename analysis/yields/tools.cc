#include "tools.h"
#include "TMath.h"

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


float minDR(const LorentzVector & lep , const vector<LorentzVector> & jet)
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

float PtMaxEta(const vector<LorentzVector> & jet)
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

std::vector<std::pair<float, std::pair<int,int> > > MOSSF(vector<Lepton> lep)
{
    int size = (int)lep.size();
    float mossf = 999999;
    std::pair<int,int> idxs = std::make_pair(-1,-1);
    std::vector<std::pair<float,std::pair<int,int>>> ret;
    if(size>=2)for(int i=0;i<size;i++){
        for(int j=i+1;j<size;j++){	
            if (is_ossf(lep[i],lep[j])) {
                float mass = (lep[i].v+lep[j].v).M();	  
                if(TMath::Abs(mass-91.2)<15.)
                {
                    mossf = mass;
                    idxs = std::make_pair(i,j);
                    ret.push_back(std::make_pair(mossf,idxs));
                }
            }	  
        }
    }

    return ret;
}

/*
 ML signal, off-Z == 8
 ML signal, on-Z == 7 
 ML fake sideband, off-Z, one L!T == 6 
 ML fake sideband, on-Z, one L!T == 5
 OS == 4
 SS signal == 3
 SS fake sideband, one L!T == 2
 SS fake sideband, two L!T == 1
*/
std::pair<int,int> get_hyp_class(std::vector<Lepton> vleps, int ss_hyp_type) {
    //
    // first let's see if we have a ML event
    // priority: off-Z, on-Z, fake sideband
    int ngoodleps = 0;
    int nfos = 0;
    for (int i = 0; i < vleps.size(); i++) {
        Lepton lep = vleps.at(i);
        if (lep.isgood and lep.ccpt > 20.) ++ngoodleps;
        if (lep.isfo and !lep.isgood and lep.ccpt > 20. and i<3) ++nfos;
    }
    std::vector<std::pair<float,std::pair<int,int>>> zcands = MOSSF(vleps); 
    int nzcands = zcands.size();

    if (ngoodleps > 0 and (ngoodleps+nfos) > 2)  {
        if (nzcands > 0) {
            if (nfos == 0) return std::make_pair(7,nzcands);
            else return std::make_pair(5,nzcands);
        }
        else {
            if (nfos == 0) return std::make_pair(8,nzcands);
            else return std::make_pair(6,nzcands);
        }
    }
    else {
        if (ss_hyp_type == 6) return std::make_pair(3,nzcands);
        else return std::make_pair(ss_hyp_type,nzcands);
    }

    return std::make_pair(-1,0);
}

bool is_ossf(Lepton lep1, Lepton lep2) {
    int idprod = lep1.id*lep2.id;
    int ret = (idprod == -121 or idprod == -169);
    return ret;
}

bool sortByPt (float i, float j) { return fabs(i-91.2) < fabs(j-91.2); }

int printLeptons(std::vector<Lepton> vleps) {
    int i = 1;
    for (auto lep : vleps) {
        std::cout << lep.name << ": pt = " << lep.v.pt() << ", eta = " << lep.v.eta() << ", is_good = " << lep.isgood << ", isfo = " << lep.isfo << std::endl;
        i++;
    }
}

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
int signal_region(int hyp_class, int nleps_tight, int nleps_fo, int nz, int nb, int nj) {
    int ret = -1;
    if (nleps_tight >=3) {
        if (nb == 1 and nj == 1) ret = 7;
        else if (nb == 1 and nj == 2) ret = 8;
        else if (nb == 1 and nj == 3) ret = 9;
        else if (nb == 1 and nj >= 4) ret = 10;
        else if (nb >= 2 and nj == 1) ret = 11;
        else if (nb >= 2 and nj == 2) ret = 12;
        else if (nb >= 2 and nj == 3) ret = 13;
        else if (nb >= 2 and nj >= 4) ret = 14;
    }
    else {
        if (nb == 1 and nj == 2) ret = 1;
        if (nb == 1 and nj == 3) ret = 2;
        if (nb == 1 and nj >= 4) ret = 3;
        if (nb >= 2 and nj == 2) ret = 4;
        if (nb >= 2 and nj == 3) ret = 5;
        if (nb >= 2 and nj >= 4) ret = 6;
    }

    return ret;
}

std::vector<std::pair<int,int>> num_z_cands(std::vector<LorentzVector> p4s, std::vector<int> ids) {
    std::vector<std::pair<int,int>> zcand_leps; 
    for (int i = 0; i < ids.size(); i++) {
        for (int j = i+1; j < ids.size(); j++) {
            if (abs(ids.at(i)) != abs(ids.at(j))) continue; // same flavor
            if (ids.at(i)*ids.at(j) > 0) continue; // OS
            LorentzVector zp4 = p4s.at(i) + p4s.at(j);
            if (fabs(zp4.mass() - 91.2) > 15) continue; // Z mass
            zcand_leps.push_back(std::make_pair(i,j));
        }
    }
    
    return zcand_leps;
}
