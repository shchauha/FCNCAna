
float electronFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.390898;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.223558;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.257473;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.275743;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0811402;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.118324;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0526458;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0393154;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0685874;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0290506;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0667149;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0931784;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0992832;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.112107;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.137475;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.10728;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0306516;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.260573;
   return 0.;
}
float electronFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.131746;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.110196;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0944638;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0605668;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0315547;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0304536;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00949893;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00569681;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00754257;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00765758;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0104451;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00912461;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0224638;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0244915;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0180825;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0635851;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0454279;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.037346;
   return 0.;
}
float electronAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.380462;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.215372;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.249431;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.269513;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0743139;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.114112;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0394612;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0311013;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0626838;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0320585;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0739721;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.036302;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0541458;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.106535;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.212249;
   return 0.;
}
float electronQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.423155;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.565624;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.406442;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.104668;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.147856;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.120556;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0372421;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0413517;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0844387;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0343948;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0593147;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0810047;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0508861;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0924386;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0796302;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.110553;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.122151;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.224711;
   return 0.;
}
float electronQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0318356;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0486777;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0671105;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0148711;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0233112;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0195014;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00668594;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00761403;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0124204;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00636132;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0129943;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0133154;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0152808;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0265066;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0164445;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0219604;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0301088;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.028383;
   return 0.;
}
float muonFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.564522;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.468891;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.647144;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.1793;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.237978;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.200655;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.047324;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0586716;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.065121;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0436406;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0364165;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.103387;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.016187;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0570286;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.117938;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0663444;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.104047;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.117938;
   return 0.;
}
float muonFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0763244;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.127832;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.233536;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0306875;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0404411;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.092672;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00425343;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00607871;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0143705;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00691032;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00799225;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0249936;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0228523;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0199916;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0663823;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0564246;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0494164;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0663823;
   return 0.;
}
float muonAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.55675;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.452566;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.638065;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.168773;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.2293;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.183889;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0321438;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0447947;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0466347;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0528542;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   return 0.;
}
float muonQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.485216;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.548741;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.662519;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.169748;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.227633;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.273098;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0569372;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0843485;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0943943;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0445833;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0516865;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0687628;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0390491;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0505665;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0450546;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0405155;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0323285;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0701499;
   return 0.;
}
float muonQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0261309;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0347082;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0587149;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00950359;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0120106;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.025437;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00317272;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00542652;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.010161;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00329927;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00448803;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0108267;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00381505;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0077086;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.00956963;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00498523;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00533713;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0238295;
   return 0.;
}
float electronFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.321581;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.279748;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.276902;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.359608;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.108563;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.172266;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0802729;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0795911;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.119644;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00716674;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0847951;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.168389;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.209688;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.183851;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.233879;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.191903;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0190088;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.398765;
   return 0.;
}
float electronFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.138818;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.142993;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.096636;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0770031;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0515106;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0431254;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0146179;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0138948;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.014755;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0238178;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0258905;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0202859;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.051457;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0499929;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0332461;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.176345;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.108516;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0539019;
   return 0.;
}
float electronAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.324251;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.281582;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.278245;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.361007;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.110654;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.173357;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0864852;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0832744;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.121746;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0296565;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.102594;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.177252;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.231916;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.20324;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.243484;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.268703;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0616552;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.40952;
   return 0.;
}
float electronQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.469227;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.522789;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.407835;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.143496;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.201734;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.146049;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0811412;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.081194;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.150573;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0798284;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.125701;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.157778;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0940064;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.172235;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.169947;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.218731;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.185943;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.282028;
   return 0.;
}
float electronQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0355347;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0714104;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0767633;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0208021;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0327217;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0254679;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0149349;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0157533;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0225891;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0148949;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0282643;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0265911;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0304175;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0477395;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0347361;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0432533;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0491912;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0372872;
   return 0.;
}
float muonFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.557524;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.350328;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.817892;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.169989;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.225226;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.295034;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0474851;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0619445;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0733372;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0418475;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0341503;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0882309;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0142134;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0684098;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.157896;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0533669;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.172812;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.157896;
   return 0.;
}
float muonFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0767693;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.124987;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.173349;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0300682;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0409736;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.103006;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0045089;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0066539;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0172637;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00797379;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00988412;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0250095;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0285866;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0275173;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0955787;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0528893;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0818124;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0955787;
   return 0.;
}
float muonAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.556319;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.347269;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.81646;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.168375;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.223807;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.293322;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0450386;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0597226;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0701198;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.033686;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0256557;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0807612;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0486041;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.12274;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00434015;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.123785;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   return 0.;
}
float muonQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.485363;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.548621;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.662505;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.169757;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.227283;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.273878;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0570671;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0852324;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0950307;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0449928;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0530749;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0711829;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0398129;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0531309;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.046544;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0419674;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0345744;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0767373;
   return 0.;
}
float muonQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0261582;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0347161;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.058717;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00950408;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.012058;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.025506;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00317994;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00549048;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0102704;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00333191;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00461377;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.011193;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00390688;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00808555;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.00991002;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00517805;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00571207;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0259197;
   return 0.;
}


float fakeRate(int id, float pt, float eta, float ht, bool isLL) { 
    if (abs(id)==11) return (isLL ? electronFakeRate(pt,eta) : electronFakeRate_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonFakeRate(pt,eta) : muonFakeRate_IsoTrigs(pt,eta));
    else return 0.;
}

float fakeRateError(int id, float pt, float eta, float ht, bool isLL) { 
    if (abs(id)==11) return (isLL ? electronFakeRateError(pt,eta) : electronFakeRateError_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonFakeRateError(pt,eta) : muonFakeRateError_IsoTrigs(pt,eta));
    else return 0.;
}

float alternativeFakeRate(int id, float pt, float eta, float ht, bool isLL) { 
    if (abs(id)==11) return (isLL ? electronAlternativeFakeRate(pt,eta) : electronAlternativeFakeRate_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonAlternativeFakeRate(pt,eta) : muonAlternativeFakeRate_IsoTrigs(pt,eta));
    else return 0.;
}

float qcdMCFakeRate(int id, float pt, float eta, float ht, bool isLL) { 
    if (abs(id)==11) return (isLL ? electronQCDMCFakeRate(pt,eta) : electronQCDMCFakeRate_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonQCDMCFakeRate(pt,eta) : muonQCDMCFakeRate_IsoTrigs(pt,eta));
    else return 0.;
}

float qcdMCFakeRateError(int id, float pt, float eta, float ht, bool isLL) { 
    if (abs(id)==11) return (isLL ? electronQCDMCFakeRateError(pt,eta) : electronQCDMCFakeRateError_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonQCDMCFakeRateError(pt,eta) : muonQCDMCFakeRateError_IsoTrigs(pt,eta));
    else return 0.;
}

