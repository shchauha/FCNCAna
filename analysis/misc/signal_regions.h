#ifndef SIGNALREGIONS_H
#define SIGNALREGIONS_H

// need to either inline everything because function bodies are in this header file
// or need to make sure it gets included only once in total
// header guard at the top only prevents against multiple inclusions in one file
// so otherwise we get the "multiple definitions of..." when compiling/linking :(
// https://stackoverflow.com/questions/8201944/multiple-definition-and-header-only-libraries

// SS stuff

enum region_t { HighHigh = 0, HighLow = 1, LowLow = 2, Multilepton = 3, LowMet = 4, Undefined = -1 };

float ptcut(int id) {
  return (abs(id) == 11) ? 15 : 10; 
}

bool is_in_lowmet(int lep1_pt, int lep2_pt, float ht, float met) {
    if (lep1_pt < 25.) return false;
    if (lep2_pt < 25.) return false;
    if (met > 50.) return false;
    if (ht < 300.) return false;
    return true;
}

region_t analysis_category_ss(int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, float ht, float met) {
    if ((lep1pt < ptcut(id1)) || (lep2pt < ptcut(id2))) return Undefined;
    // Want 25, 15/10, 15/10 (babymaker guarantees everything is above 15/10)
    // meaning >=1 of the 3 has to be >25 and >=2 have to be >20
    int numhigh = (lep1pt > 25.) + (lep2pt > 25.);

    // New categories
    if ( nleps > 2
            and (numhigh>=1)
            // and ((lep1pt>25.)+(lep2pt>25.)+(lep3pt>25.))>=1
            // and ((lep1pt>10.)+(lep2pt>10.)+(lep3pt>10.))>=2
       ) return Multilepton;
    if (numhigh == 2) {
        if (is_in_lowmet(lep1pt,lep2pt,ht,met)) return LowMet;
        else return HighHigh;
    }
    if (numhigh == 1) return HighLow;
    return LowLow;

    // // 2016 categories
    // if (numhigh == 2) return HighHigh;
    // if (numhigh == 1) return HighLow;
    // return LowLow;

}


bool passes_baseline_ss(int njets, int nbtags, float met, float ht, int id1, int id2, float lep1_pt, float lep2_pt) {
  //Kinematic cuts
  if (lep1_pt < ptcut(id1)) return false;
  if (lep2_pt < ptcut(id2)) return false;
  if (njets < 2) return false;
  if ((met < 50.) and not is_in_lowmet(lep1_pt,lep2_pt,ht,met)) return false;
  return true;
}

int getNsrsLM() { return 11; }
int signal_region_lowmet(int njets, int nbtags, float met, float ht, float lep1pt, float lep2pt) {
    // pg 12 of https://indico.cern.ch/event/766192/contributions/3180522/attachments/1735903/2807690/20181019_RA5_LucienLo_draft.pdf
    if (!is_in_lowmet(lep1pt, lep2pt, ht, met)) return -1;
    if (ht < 1125.) {
        if (nbtags == 0) {
            if (njets < 5) return 1;
            else return 2;
        } else if (nbtags == 1) {
            if (njets < 5) return 3;
            else return 4;
        } else if (nbtags == 2) {
            if (njets < 5) return 5;
            else return 6;
        } else { // pg 32 of 2016AN. only one njet bin for nb==3 http://uaf-8.t2.ucsd.edu/~namin/dump/AN-16-386_temp.pdf
            return 7;
        }
    }
    else if (ht < 1300.) {
        if (njets < 5) return 8;
        else return 9;
    }
    else {
        if (njets < 5) return 10;
        else return 11;
    }
    return -1;
}

int getNsrsML() { return 44; }
int signal_region_multilepton(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, bool isClass6) {
    if (met < 50) return -1; 
    // pg 6 of https://arxiv.org/pdf/1710.09154.pdf
    // note mtmin is the mt of the lepton that doesnt' make the z if it's a class6 event
    // off z is 1-23; onz is 24-46
    int imt = mt_min > 120.;
    bool onz = isClass6;
    if (not onz) {
        if (met > 300) return 20+imt;
        if (ht > 600) {
            if (met > 150) return 18+imt;
            else return 16+imt;
        }
        if (nbtags >= 3) return 15;
        if (nbtags == 0 and ht < 400 and met < 150) return 1+imt;
        if (nbtags == 0 and ht < 400 and met > 150) return 3+imt;
        if (nbtags == 0 and ht > 400 and met < 150) return 5;
        if (nbtags == 0 and ht > 400 and met > 150) return 6;
        if (nbtags == 1 and ht < 400 and met < 150) return 7;
        if (nbtags == 1 and ht < 400 and met > 150) return 8;
        if (nbtags == 1 and ht > 400 and met < 150) return 9;
        if (nbtags == 1 and ht > 400 and met > 150) return 10;
        if (nbtags == 2 and ht < 400 and met < 150) return 11;
        if (nbtags == 2 and ht < 400 and met > 150) return 12;
        if (nbtags == 2 and ht > 400 and met < 150) return 13;
        if (nbtags == 2 and ht > 400 and met > 150) return 14;
    } else {
        if (met > 300) return 43+imt;
        if (ht > 600) {
            if (met > 150) return 41+imt;
            else return 39+imt;
        }
        if (nbtags >= 3) return 38;
        if (nbtags == 0 and ht < 400 and met < 150) return 22+imt;
        if (nbtags == 0 and ht < 400 and met > 150) return 24+imt;
        if (nbtags == 0 and ht > 400 and met < 150) return 26+imt;
        if (nbtags == 0 and ht > 400 and met > 150) return 28+imt;
        if (nbtags == 1 and ht < 400 and met < 150) return 30;
        if (nbtags == 1 and ht < 400 and met > 150) return 31;
        if (nbtags == 1 and ht > 400 and met < 150) return 32;
        if (nbtags == 1 and ht > 400 and met > 150) return 33;
        if (nbtags == 2 and ht < 400 and met < 150) return 34;
        if (nbtags == 2 and ht < 400 and met > 150) return 35;
        if (nbtags == 2 and ht > 400 and met < 150) return 36;
        if (nbtags == 2 and ht > 400 and met > 150) return 37;
    }
    return -1;
}

int getNsrsHH() { return 62; }
int getNsrsHL() { return 43; }
int getNsrsLL() { return 8; }
int signal_region_ss(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, bool isClass6, float mtnonz) {

    //Calculate lep_pt
    region_t lep_pt = analysis_category_ss(id1, id2, lep1pt, lep2pt, lep3pt, nleps, ht, met); 

    if (lep_pt == Multilepton) {
        // Note mtnonz instead of mt_min
        return signal_region_multilepton(njets, nbtags, met, ht, mtnonz, id1, id2, lep1pt, lep2pt, lep3pt, nleps, isClass6);
    }

    if (lep_pt == LowMet and not isClass6 and nleps == 2) {
        // according to https://twiki.cern.ch/twiki/bin/view/CMS/RA5Run2
        // no low met region if the event is onZ, or even if it just has a 3rd tight lepton
        int lmsr = signal_region_lowmet(njets, nbtags, met, ht, lep1pt, lep2pt);
        if (lmsr > 0) return lmsr;
    }

    // After low met region
    if (met < 50) return -1; 

    // remember that sgn(pdgid) != sgn(charge), it's flipped. so mad.
    int mm = (id1 > 0);

    if (lep_pt != LowLow && met>300 && ht<300) return -1;

    //High-high
    if (lep_pt == HighHigh) {
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
	  if (nbtags == 0){ 
	    //cout<< "nb "<<nbtags<<" mt "<<mt_min<<" met "<<met<<" njets "<<njets<<" ht "<<ht<<" l1pt "<<lep1pt<<" l2pt "<<lep2pt<<endl;	      
	    return 3;
	  } 
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
    }

    // High-Low
    if (lep_pt == HighLow){
        if (met >= 300 && ht >= 300) {
            if (njets <= 4) {
                if (met < 500) return 33+mm;
                else return 35+mm;
            } else {
                if (met < 500) return 37+mm;
                else return 39;
            }
        }
        if (ht >= 1125) {
            if (ht < 1300) return 40+mm;
            else return 42+mm;
        }
        if (ht < 300){ 
            if (mt_min >= 120) return 31;
            if (nbtags == 0 && met < 200 && njets <= 4) return 1; 
            if (nbtags == 0) return 3;
            if (nbtags == 1 && met < 200 && njets <= 4) return 8; 
            if (nbtags == 1) return 10+mm;
            if (nbtags == 2 && met < 200 && njets <= 4) return 17; 
            if (nbtags == 2) return 19+mm;
            if (nbtags >= 3) return 26+mm;
        }  
        if (ht >= 300){
            if (mt_min >= 120) return 32;
            if (nbtags == 0 && met < 200 && njets <= 4) return 2; 
            if (nbtags == 0 && met < 200 && njets > 4) return 4; 
            if (nbtags == 0 && met < 300 && njets <= 4) return 5+mm; 
            if (nbtags == 0 && met < 300 && njets > 4) return 7; 
            if (nbtags == 1 && met < 200 && njets <= 4) return 9; 
            if (nbtags == 1 && met < 200 && njets > 4) return 12+mm; 
            if (nbtags == 1 && met < 300 && njets <= 4) return 14;
            if (nbtags == 1 && met < 300 && njets > 4) return 15+mm; 
            if (nbtags == 2 && met < 200 && njets <= 4) return 18; 
            if (nbtags == 2 && met < 200 && njets > 4) return 21+mm; 
            if (nbtags == 2 && met < 300 && njets <= 4) return 23+mm; 
            if (nbtags == 2 && met < 300 && njets > 4) return 25; 
            if (nbtags >= 3 && met < 200) return 28+mm; 
            if (nbtags >= 3 && met >= 200) return 30;
        }
    }

  // Low-Low
  if (lep_pt == LowLow) {
      if (ht < 400) return -1;  // Raise threshold since turn on worse in 2017/2018
      if (mt_min > 120) return 8; 
      if (nbtags == 0 && met < 200) return 1;
      if (nbtags == 0 && met >= 200) return 2;
      if (nbtags == 1 && met < 200) return 3;
      if (nbtags == 1 && met >= 200) return 4;
      if (nbtags == 2 && met < 200) return 5;
      if (nbtags == 2 && met >= 200) return 6;
      if (nbtags >= 3) return 7;
  }

  // Otherwise undefined
  cout << "WARNING: SR UNDEFINED (should never get here)" << endl;
  cout << "  --> lepton pts are: " << lep1pt << " " << lep2pt << endl;
  cout << "  --> ht & met are: " << ht << " " << met << endl;
  cout << "  --> njets & nbtags: " << njets << " " << nbtags << endl;
  return -1;
}


bool passes_baseline_ft(int njets, int nbtags, float met, float ht, int id1, int id2, float lep1_pt, float lep2_pt) {
    if (lep1_pt < 25.) return 0;
    else if (lep2_pt < 20.) return 0;
    else if (njets < 2) return 0;
    else if (nbtags < 2) return 0;
    else if (ht < 300) return 0;
    else if (met < 50) return 0;
    else return 1;
}

bool passes_baseline_ss(int njets, int nbtags, float met, float ht, int id1, int id2, float lep1_pt, float lep2_pt, float metcut) {
  if (lep1_pt < ptcut(id1)) return false;
  if (lep2_pt < ptcut(id2)) return false;
  if (njets < 2) return false;
  if (met < metcut) return false;
  return true;
}


/*
# paste the below into the python interpreter
import numpy as np
def get_sigmoid_binning(nbins):
    x = np.linspace(0.,1.,nbins+1)
    # k = 11.20-0.22*nbins # XXX rederive every time we re-train
    k = 11.5 # test effect of flat k independent of nbins
    k = k*np.ones(len(x))
    k[x < 0.5] *= 0.5 # for bkg enriched, don't split too finely, so soften the sigmoid by 1/2 -- threshold of 0.5 to keep continuity
    bins = 1.0/(1+np.exp(-k*(x-0.5)))+2.0*(x-0.5)/(1+np.exp(k/2))
    return bins

# note that the number of SRs is actually 1 higher than `nbins` because we add CRZ in front
nbins=13; print "int getBDTBin(float disc, bool crz) {\n    if (crz) return 1;\n"+"\n".join(["    if (disc > {:.4f}) return {};".format(x,ix+2) for ix,x in list(enumerate(get_sigmoid_binning(nbins)[:-1]))[::-1]])+"\n    return -1;\n}}\nint getNsrsDisc() {{ return {}; }}".format(nbins+1)
# or for many at once, labeled by the bin numbers
for nbins in range(10,21): print "int getBDTBin_{}(float disc, bool crz) {{\n    if (crz) return 1;\n".format(nbins)+"\n".join(["    if (disc > {:.4f}) return {};".format(x,ix+2) for ix,x in list(enumerate(get_sigmoid_binning(nbins)[:-1]))[::-1]])+"\n    return -1;\n}}\nint getNsrsDisc_{}() {{ return {}; }}".format(nbins,nbins+1)
*/

int getBDTBin_17(float disc, bool crz) {
    if (crz) return 1;

    ////old
    //if (disc > 0.9966) return 18;
    //if (disc > 0.9903) return 17;
    //if (disc > 0.9784) return 16;
    //if (disc > 0.9562) return 15;
    //if (disc > 0.9156) return 14;
    //if (disc > 0.8453) return 13;
    //if (disc > 0.7345) return 12;
    //if (disc > 0.5839) return 11;
    //if (disc > 0.4547) return 10;
    //if (disc > 0.3664) return 9;
    //if (disc > 0.2847) return 8;
    //if (disc > 0.2124) return 7;
    //if (disc > 0.1509) return 6;
    //if (disc > 0.1001) return 5;
    //if (disc > 0.0590) return 4;
    //if (disc > 0.0262) return 3;
    //if (disc > 0.0000) return 2;

    if (disc > 0.9956) return 18;
    if (disc > 0.9884) return 17;
    if (disc > 0.9749) return 16;
    if (disc > 0.9506) return 15;
    if (disc > 0.9076) return 14;
    if (disc > 0.8357) return 13;
    if (disc > 0.7260) return 12;
    if (disc > 0.6054) return 11;
    if (disc > 0.4741) return 10;
    if (disc > 0.3704) return 9;
    if (disc > 0.2905) return 8;
    if (disc > 0.2190) return 7;
    if (disc > 0.1573) return 6;
    if (disc > 0.1055) return 5;
    if (disc > 0.0659) return 4;
    if (disc > 0.0362) return 3;
    if (disc > 0.0000) return 2;

    return -1;
}
int getNsrsDisc_17() { return 18; }

int getNsrsDisc(int nbins) {
    if (nbins == 17) return getNsrsDisc_17();
    return -1;
}
int getBDTBin(int nbins, float disc, bool crz) {
    if (nbins == 17) return getBDTBin_17(disc,crz);
    return -1;
}


// // nominal
// int signal_region_ft(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, int isClass6){
//     if (lep1pt < 25.) return -1;
//     if (lep2pt < 20.) return -1;
//     if (njets < 2) return -1;
//     if (nbtags < 2) return -1;
//     if (ht < 300) return -1;
//     if (met < 50) return -1;
//     if (isClass6) {
//         if (lep3pt >= 20.) return 1;
//         else return -1;
//     }
//     if (nleps == 2) {
//         if (nbtags == 2) {
//             if (njets <= 5) return 2;
//             if (njets == 6) return 3;
//             if (njets == 7) return 4;
//             if (njets >= 8) return 5;
//         } else if (nbtags == 3) {
//             if (njets == 5) return 6;
//             if (njets == 6) return 6;
//             if (njets >= 7) return 7;
//         } else if (nbtags >= 4) {
//             if (njets >= 5) return 8;
//         }
//     } else {
//         if (nbtags == 2) {
//             if (njets >= 5) return 9;
//         } else if (nbtags >= 3) {
//             if (njets >= 4) return 10;
//         }
//     }
//     return -1;
// }
// int getNsrs() { return 10; } // note that this includes CRs
// int getNCRs() { return 2; }

// // 18 bins
// int signal_region_ft(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, bool isClass6){
//     if (lep1pt < 25.) return -1;
//     if (lep2pt < 20.) return -1;
//     if (njets < 2) return -1;
//     if (nbtags < 2) return -1;
//     if (ht < 300) return -1;
//     if (met < 50) return -1;
//     if (isClass6) {
//         if (lep3pt >= 20.) return 1;
//         else return -1;
//     }
//     if (nleps == 2) {
//         if (nbtags == 2) {
//             if (njets <= 5) return 2;
//             if (njets == 6) return 3;
//             if (njets == 7) return 4;
//             if (njets >= 8) return 5;
//         } else if (nbtags == 3) {
//             if (njets == 5) return 6;
//             if (njets == 6) return 7;
//             if (njets == 7) return 8;
//             if (njets >= 8) return 9;
//         } else if (nbtags >= 4) {
//             if (njets == 5) return 10;
//             if (njets == 6) return 11;
//             if (njets >= 7) return 12;
//         }
//     } else {
//         if (nbtags == 2) {
//             if (njets == 5) return 13;
//             if (njets == 6) return 14;
//             if (njets >= 7) return 15;
//         } else if (nbtags >= 3) {
//             if (njets == 4) return 16;
//             if (njets == 5) return 17;
//             if (njets >= 6) return 18;
//         }
//     }
//     return -1;
// }
// int getNsrs() { return 18; } // note that this includes CRs
// int getNCRs() { return 2; }
// // int getNsrsDisc() { return 13; }

// 16 bins -- MERGING the three 4-b fake dominated bins
int signal_region_ft(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, bool isClass6){
    if (lep1pt < 25.) return -1;
    if (lep2pt < 20.) return -1;
    if (njets < 2) return -1;
    if (nbtags < 2) return -1;
    if (ht < 300) return -1;
    if (met < 50) return -1;
    if (isClass6) {
        if (lep3pt >= 20.) return 1;
        else return -1;
    }
    if (nleps == 2) {
        if (nbtags == 2) {
            if (njets <= 5) return 2;
            if (njets == 6) return 3;
            if (njets == 7) return 4;
            if (njets >= 8) return 5;
        } else if (nbtags == 3) {
            if (njets == 5) return 6;
            if (njets == 6) return 7;
            if (njets == 7) return 8;
            if (njets >= 8) return 9;
        } else if (nbtags >= 4) {
            if (njets >= 5) return 10;
        }
    } else {
        if (nbtags == 2) {
            if (njets == 5) return 11;
            if (njets == 6) return 12;
            if (njets >= 7) return 13;
        } else if (nbtags >= 3) {
            if (njets == 4) return 14;
            if (njets == 5) return 15;
            if (njets >= 6) return 16;
        }
    }
    return -1;
}
int getNsrs() { return 16; } // note that this includes CRs
int getNCRs() { return 2; }
// int getNsrsDisc() { return 13; }

// // 16 bins -- DROPPING the three 4-b fake dominated bins
// int signal_region_ft(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, bool isClass6){
//     if (lep1pt < 25.) return -1;
//     if (lep2pt < 20.) return -1;
//     if (njets < 2) return -1;
//     if (nbtags < 2) return -1;
//     if (ht < 300) return -1;
//     if (met < 50) return -1;
//     if (isClass6) {
//         if (lep3pt >= 20.) return 1;
//         else return -1;
//     }
//     if (nleps == 2) {
//         if (nbtags == 2) {
//             if (njets <= 5) return 2;
//             if (njets == 6) return 3;
//             if (njets == 7) return 4;
//             if (njets >= 8) return 5;
//         } else if (nbtags == 3) {
//             if (njets == 5) return 6;
//             if (njets == 6) return 7;
//             if (njets == 7) return 8;
//             if (njets >= 8) return 9;
//         } else if (nbtags >= 4) {
//             return -1;
//         }
//     } else {
//         if (nbtags == 2) {
//             if (njets == 5) return 11;
//             if (njets == 6) return 12;
//             if (njets >= 7) return 13;
//         } else if (nbtags >= 3) {
//             if (njets == 4) return 14;
//             if (njets == 5) return 15;
//             if (njets >= 6) return 16;
//         }
//     }
//     return -1;
// }
// int getNsrs() { return 16; } // note that this includes CRs
// int getNCRs() { return 2; }
// // int getNsrsDisc() { return 13; }


// // 19 bins  with q split CRW
// int signal_region_ft(int njets, int nbtags, float met, float ht, float mt_min, int id1, int id2, float lep1pt, float lep2pt, float lep3pt, int nleps, int isClass6){
//     if (lep1pt < 25.) return -1;
//     if (lep2pt < 20.) return -1;
//     if (njets < 2) return -1;
//     if (nbtags < 2) return -1;
//     if (ht < 300) return -1;
//     if (met < 50) return -1;
//     if (isClass6) {
//         if (lep3pt >= 20.) return 1;
//         else return -1;
//     }
//     if (nleps == 2) {
//         if (nbtags == 2) {
//             if (njets <= 5) return 2 + (id1>0);
//             if (njets == 6) return 4;
//             if (njets == 7) return 5;
//             if (njets >= 8) return 6;
//         } else if (nbtags == 3) {
//             if (njets == 5) return 7;
//             if (njets == 6) return 8;
//             if (njets == 7) return 9;
//             if (njets >= 8) return 10;
//         } else if (nbtags >= 4) {
//             if (njets == 5) return 11;
//             if (njets == 6) return 12;
//             if (njets >= 7) return 13;
//         }
//     } else {
//         if (nbtags == 2) {
//             if (njets == 5) return 14;
//             if (njets == 6) return 15;
//             if (njets >= 7) return 16;
//         } else if (nbtags >= 3) {
//             if (njets == 4) return 17;
//             if (njets == 5) return 18;
//             if (njets >= 6) return 19;
//         }
//     }
//     return -1;
// }
// int getNsrs() { return 19; } // note that this includes CRs
// int getNCRs() { return 3; }
// int getNsrsDisc() { return 13; }


#endif
