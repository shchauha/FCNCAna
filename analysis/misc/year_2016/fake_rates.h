float electronFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.212947;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.191067;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.114535;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0598557;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0583876;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0485127;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0502397;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.055552;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0547698;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0612282;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0709114;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0685393;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.126793;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.134804;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.138614;
   return 0.;
}
float electronFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00957134;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00844421;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00823802;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00347422;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00319052;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00318651;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00200887;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0019377;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00178593;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00440903;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00399083;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00302364;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0125623;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0100297;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00661812;
   return 0.;
}
float electronAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.213998;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.191629;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.115074;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0615844;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0593789;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0493942;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0560121;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0587569;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0571373;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0849149;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0881262;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0790308;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.178784;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.167336;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.155869;
   return 0.;
}
float electronQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.300186;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.294154;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.249293;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.08778;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0977887;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0902851;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0540886;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0877336;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0624113;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0590063;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.111845;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0733026;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.145939;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.161412;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.14486;
   return 0.;
}
float muonFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.303083;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.358831;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.383154;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0711473;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0771226;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0880371;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0423236;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0447963;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0513743;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0452653;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0421763;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0437186;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0587507;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0498613;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   return 0.;
}
float muonFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00721329;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00967646;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0206302;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00302625;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00393619;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.00914311;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00091705;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00116013;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.00304744;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00166003;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00196206;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0048238;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00547084;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0062721;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0136565;
   return 0.;
}
float muonAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.30331;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.359038;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.383389;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0715748;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0775189;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0884769;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0435583;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0458865;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0529224;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0498163;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0467553;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0499567;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.081462;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0718956;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0185418;
   return 0.;
}
float muonQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.351064;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.400847;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.442068;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0909698;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.101132;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.127172;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0562548;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0566511;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0691028;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0554049;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.048622;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0559193;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.056687;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0457155;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0481817;
   return 0.;
}

float electronFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.29519;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.259057;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.131562;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.12745;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.127191;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0859382;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.14444;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.146665;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.129899;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.197826;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.196305;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.171395;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.328461;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.298701;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.271449;
   return 0.;
}
float electronFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0118739;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0112231;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00964814;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00706781;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00665043;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00555776;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00618466;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00555629;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0044183;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0126553;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0104359;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00714557;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0283268;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0210151;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0120671;
   return 0.;
}
float electronAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.29629;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.259744;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.132181;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.130749;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.129121;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0874116;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.159532;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.154472;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.13482;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.256636;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.236685;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.193891;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.428686;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.359844;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.300518;
   return 0.;
}
float electronQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.376709;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.368919;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.257943;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.16129;
   if (pt>=15 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.174033;
   if (pt>=15 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.149289;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.129813;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.162479;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.103529;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.134471;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.196851;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.109899;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.284008;
   if (pt>=50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.277981;
   if (pt>=50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.220997;
   return 0.;
}
float muonFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.334722;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.372321;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.377708;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0913886;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.101491;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.106955;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0620827;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0686233;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0652329;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0691848;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0710012;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0571331;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.087843;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.089634;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   return 0.;
}
float muonFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00545534;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00709176;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0157604;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00273855;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00363839;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.00773028;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00128476;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00173027;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.00356143;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00254016;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00329751;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.00663705;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00883208;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0113382;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0306071;
   return 0.;
}
float muonAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.334989;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.372567;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.378029;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0920142;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.102081;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.107549;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0641541;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0705977;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0674937;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0772773;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0797381;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0671027;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.128024;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.132056;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0410995;
   return 0.;
}

float muonQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.370786;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.421886;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.455732;
   if (pt>=15 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.115238;
   if (pt>=15 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.126395;
   if (pt>=15 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.149923;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0838494;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0826661;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0929938;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0852247;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0755065;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0774966;
   if (pt>=50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.093793;
   if (pt>=50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0743385;
   if (pt>=50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0753781;
   return 0.;
}

float electronFakeRate_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2331; // +-6.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2390; // +-7.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2797; // +-6.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0821; // +-17.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0726; // +-11.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1338; // +-8.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0482; // +-6.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0624; // +-5.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0952; // +-4.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0543; // +-13.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0560; // +-11.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0954; // +-6.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0715; // +-28.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0945; // +-17.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1443; // +-10.4%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1368; // +-41.8%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1362; // +-35.3%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1962; // +-17.5%;
   return 0.;
}
float electronFakeRateError_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0156;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0167;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0186;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0144;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0086;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0114;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0033;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0034;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0041;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0075;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0063;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0061;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0203;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0169;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0151;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0572;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0481;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0343;
   return 0.;
}

float electronAlternativeFakeRate_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2362; // +-6.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2414; // +-6.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2814; // +-6.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0859; // +-16.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0752; // +-11.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1359; // +-8.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0565; // +-5.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0674; // +-5.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0991; // +-4.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0861; // +-8.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0785; // +-7.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1116; // +-5.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1276; // +-14.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1309; // +-12.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1704; // +-8.5%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2483; // +-18.7%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2220; // +-18.3%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2458; // +-12.8%;
   return 0.;
}
float electronAlternativeFakeRateError_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0155;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0167;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0185;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0143;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0086;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0113;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0033;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0034;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0041;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0072;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0061;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0061;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0185;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0160;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0145;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0465;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0407;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0314;
   return 0.;
}

float electronQCDMCFakeRate_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2995; // +-3.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3090; // +-4.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3295; // +-6.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0906; // +-7.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1137; // +-10.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1105; // +-10.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0482; // +-8.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0612; // +-9.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0776; // +-10.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0631; // +-14.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0771; // +-15.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0820; // +-11.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1278; // +-18.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1032; // +-19.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1244; // +-15.3%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1306; // +-25.1%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1174; // +-13.0%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2522; // +-13.9%;
   return 0.;
}
float electronQCDMCFakeRateError_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0091;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0128;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0210;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0065;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0114;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0111;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0040;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0059;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0079;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0094;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0117;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0090;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0238;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0197;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0190;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0328;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0152;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0350;
   return 0.;
}

float muonFakeRate_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3401; // +-3.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3969; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.4729; // +-6.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0996; // +-7.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1058; // +-8.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1425; // +-16.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0422; // +-13.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0645; // +-13.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1249; // +-20.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0410; // +-3.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0458; // +-4.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0762; // +-8.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0427; // +-6.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0426; // +-8.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0684; // +-17.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0431; // +-20.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0540; // +-20.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0360; // +-66.0%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1002; // +-30.6%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0701; // +-54.7%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0360; // +-66.0%;
   return 0.;
}
float muonFakeRateError_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0110;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0145;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0299;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0069;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0085;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0232;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0058;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0087;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0259;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0015;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0020;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0068;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0029;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0034;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0118;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0089;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0111;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0237;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0307;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0383;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0237;
   return 0.;
}

float muonAlternativeFakeRate_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3403; // +-3.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3970; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.4730; // +-6.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0999; // +-6.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1061; // +-8.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1428; // +-16.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0427; // +-13.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0650; // +-13.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1255; // +-20.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0422; // +-3.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0470; // +-4.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0779; // +-8.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0475; // +-6.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0475; // +-7.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0763; // +-15.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0610; // +-14.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0709; // +-15.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0554; // +-44.4%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1503; // +-20.6%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1280; // +-28.7%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0378; // +-149.4%;
   return 0.;
}
float muonAlternativeFakeRateError_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0110;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0145;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0299;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0069;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0085;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0232;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0058;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0087;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0259;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0015;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0020;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0069;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0029;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0035;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0122;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0091;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0112;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0246;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0310;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0368;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0564;
   return 0.;
}

float muonQCDMCFakeRate_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3730; // +-1.1%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.4086; // +-1.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.4618; // +-2.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1071; // +-1.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1230; // +-1.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1691; // +-2.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0643; // +-2.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0725; // +-2.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1068; // +-4.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0535; // +-2.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0595; // +-2.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0803; // +-4.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0507; // +-4.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0505; // +-4.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0704; // +-8.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0483; // +-8.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0421; // +-10.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0681; // +-17.6%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0423; // +-17.9%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0495; // +-19.2%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0483; // +-40.2%;
   return 0.;
}
float muonQCDMCFakeRateError_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0040;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0053;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0098;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0015;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0020;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0047;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0013;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0018;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0042;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0012;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0016;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0037;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0020;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0024;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0058;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0042;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0046;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0120;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0076;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0095;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0194;
   return 0.;
}

float electronFakeRate_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3017; // +-5.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2834; // +-6.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3555; // +-5.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1231; // +-9.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1661; // +-8.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1847; // +-6.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1141; // +-6.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1455; // +-5.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1769; // +-4.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1586; // +-12.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1496; // +-10.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2030; // +-5.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2071; // +-25.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2136; // +-17.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2754; // +-9.1%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3802; // +-26.1%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2779; // +-29.0%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2837; // +-15.8%;
   return 0.;
}
float electronFakeRateError_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0176;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0194;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0179;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0122;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0133;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0125;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0078;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0080;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0072;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0197;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0158;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0120;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0526;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0364;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0252;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0991;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0805;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0449;
   return 0.;
}

float electronAlternativeFakeRate_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3057; // +-5.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2871; // +-6.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3573; // +-5.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1296; // +-9.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1700; // +-7.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1874; // +-6.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1310; // +-5.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1555; // +-5.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1829; // +-3.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2268; // +-7.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1973; // +-7.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2300; // +-5.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3218; // +-12.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2792; // +-11.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3113; // +-7.5%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.5223; // +-12.4%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3946; // +-15.3%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3387; // +-11.8%;
   return 0.;
}
float electronAlternativeFakeRateError_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0175;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0193;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0178;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0121;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0133;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0124;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0077;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0079;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0072;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0172;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0144;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0115;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0405;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0318;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0234;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0650;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0604;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0400;
   return 0.;
}

float electronQCDMCFakeRate_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3789; // +-3.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3942; // +-4.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.4086; // +-6.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1643; // +-7.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1805; // +-10.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1628; // +-10.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1226; // +-8.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1298; // +-10.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1335; // +-10.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1568; // +-15.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1756; // +-14.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1381; // +-11.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2846; // +-17.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1756; // +-20.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2120; // +-15.1%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2561; // +-25.1%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2018; // +-14.6%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3416; // +-13.2%;
   return 0.;
}
float electronQCDMCFakeRateError_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0115;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0159;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0250;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0120;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0186;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0167;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0105;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0136;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0142;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0240;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0259;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0158;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0491;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0363;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0320;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0644;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0295;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0450;
   return 0.;
}

float muonFakeRate_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3816; // +-2.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.4318; // +-2.5%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.4678; // +-4.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1133; // +-5.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1290; // +-5.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1554; // +-11.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0789; // +-8.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0940; // +-9.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1229; // +-17.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0585; // +-3.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0697; // +-4.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0927; // +-8.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0668; // +-6.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0662; // +-7.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0916; // +-16.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0704; // +-21.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0829; // +-21.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0723; // +-78.5%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1255; // +-32.9%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0928; // +-63.0%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0723; // +-78.5%;
   return 0.;
}
float muonFakeRateError_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0084;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0108;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0226;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0057;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0077;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0176;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0064;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0085;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0210;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0019;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0028;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0074;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0042;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0051;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0147;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0149;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0175;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0567;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0413;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0585;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0567;
   return 0.;
}

float muonAlternativeFakeRate_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3818; // +-2.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.4320; // +-2.5%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.4680; // +-4.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1137; // +-5.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1295; // +-5.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1559; // +-11.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0797; // +-8.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0948; // +-8.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1238; // +-16.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0604; // +-3.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0717; // +-4.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0952; // +-7.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0756; // +-5.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0753; // +-6.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1044; // +-14.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1062; // +-13.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1145; // +-14.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1293; // +-39.6%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2047; // +-18.7%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1927; // +-25.8%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0625; // +-146.6%;
   return 0.;
}
float muonAlternativeFakeRateError_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0084;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0108;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0225;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0057;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0077;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0176;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0064;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0084;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0209;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0019;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0028;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0075;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0043;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0052;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0149;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0144;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0169;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0511;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0383;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0497;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0916;
   return 0.;
}

float muonQCDMCFakeRate_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3940; // +-1.1%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.4297; // +-1.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.4809; // +-2.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1301; // +-1.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1486; // +-1.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1959; // +-2.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0885; // +-2.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0979; // +-2.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1375; // +-4.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0805; // +-2.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0878; // +-2.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1090; // +-4.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0785; // +-3.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0788; // +-4.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1008; // +-8.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0793; // +-8.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0707; // +-10.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1073; // +-17.2%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0748; // +-17.7%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0813; // +-19.4%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0930; // +-39.3%;
   return 0.;
}
float muonQCDMCFakeRateError_IsoTrigs_legacy(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0042;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0055;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0102;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0025;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0055;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0024;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0054;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0024;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0050;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0031;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0037;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0083;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0068;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0077;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0184;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0132;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0158;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0365;
   return 0.;
}


float fakeRate(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
    // SS uses iso for 2016; FT still uses nonisoHT trigs for ht>300 in 2016
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronFakeRate(pt,eta);
    else if (abs(id)==13) return muonFakeRate(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronFakeRate_legacy(pt,eta) : electronFakeRate_IsoTrigs_legacy(pt,eta));
    else if (abs(id)==13) return (isLL ? muonFakeRate_legacy(pt,eta) : muonFakeRate_IsoTrigs_legacy(pt,eta));
    else return 0.;
  }
}

float fakeRateError(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronFakeRateError(pt,eta);
    else if (abs(id)==13) return muonFakeRateError(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronFakeRateError_legacy(pt,eta) : electronFakeRateError_IsoTrigs_legacy(pt,eta));
    else if (abs(id)==13) return (isLL ? muonFakeRateError_legacy(pt,eta) : muonFakeRateError_IsoTrigs_legacy(pt,eta));
    else return 0.;
  }
}

float alternativeFakeRate(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronAlternativeFakeRate(pt,eta);
    else if (abs(id)==13) return muonAlternativeFakeRate(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronAlternativeFakeRate_legacy(pt,eta) : electronAlternativeFakeRate_IsoTrigs_legacy(pt,eta));
    else if (abs(id)==13) return (isLL ? muonAlternativeFakeRate_legacy(pt,eta) : muonAlternativeFakeRate_IsoTrigs_legacy(pt,eta));
    else return 0.;
  }
}

float qcdMCFakeRate(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronQCDMCFakeRate(pt,eta);
    else if (abs(id)==13) return muonQCDMCFakeRate(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronQCDMCFakeRate_legacy(pt,eta) : electronQCDMCFakeRate_IsoTrigs_legacy(pt,eta));
    else if (abs(id)==13) return (isLL ? muonQCDMCFakeRate_legacy(pt,eta) : muonQCDMCFakeRate_IsoTrigs_legacy(pt,eta));
    else return 0.;
  }
}

float qcdMCFakeRateError(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
    return 0.01;
}

