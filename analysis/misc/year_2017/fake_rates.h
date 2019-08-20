
float electronFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.369799;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.383711;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.362114;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.132995;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0755862;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.13318;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0315842;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0375936;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0638441;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0251564;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0372046;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0549404;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0626212;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0735135;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0930926;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0615474;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0921124;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.215674;
   return 0.;
}
float electronFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.111517;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.135131;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.117339;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0446104;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.030019;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0434261;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00535749;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00487832;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00649873;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00893807;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00881641;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.00836383;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.02068;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0186961;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.015283;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0546291;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0543308;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0320888;
   return 0.;
}
float electronAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.370706;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.384237;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.362589;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.133573;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.076063;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.133558;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0323423;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0380301;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0642812;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.027257;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0390489;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0563199;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0655426;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0756125;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0946766;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0687023;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0973192;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.218059;
   return 0.;
}
float electronQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.419199;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.555941;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.401207;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.114011;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.158569;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.133304;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0314069;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0421826;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0801429;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0316121;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0632297;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0760515;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0463065;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0823888;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.089693;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.116887;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.109561;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.220749;
   return 0.;
}
float electronQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0342985;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0555685;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0727932;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0205042;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0282231;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0253392;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.00515412;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.00897237;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0128846;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.006374;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0146773;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0143184;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0170263;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0251064;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0198553;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0240659;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0264197;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0285809;
   return 0.;
}
float muonFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.354001;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.621269;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.867166;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.116507;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.158847;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.234055;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0543563;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0777358;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.106282;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0327625;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0463964;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0733981;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0309343;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0790267;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.171842;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0407916;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0670991;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0746806;
   return 0.;
}
float muonFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0698244;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0833711;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0705841;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0251053;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0322079;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0688975;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00340536;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00685065;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0130358;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0049077;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00615799;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0156267;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.012758;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0173773;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0494579;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0456848;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0493683;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0821201;
   return 0.;
}
float muonAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.35888;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.623172;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.867347;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.12063;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.163059;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.237723;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0592251;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0819246;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.110911;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0470295;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0591944;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.08813;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0646151;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.110517;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.212474;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.14143;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.148193;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.139988;
   return 0.;
}
float muonQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.48737;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.546696;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.623436;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.180422;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.233761;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.276283;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0559878;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0855424;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0870273;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0455401;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0476689;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0605103;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.039713;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0502564;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0467387;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0415329;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0309944;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0745517;
   return 0.;
}
float muonQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0285992;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0377265;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0656202;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0106153;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0135292;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.027727;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00343024;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0062243;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0105304;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00374776;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00434191;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0107859;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00430192;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00834368;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.010489;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00565895;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00597259;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0284783;
   return 0.;
}
float electronFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.399488;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.30978;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.434417;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.134325;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0698102;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.176644;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0666072;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0726379;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.127929;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0813633;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0947138;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.117984;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.200698;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.211528;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.202261;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.179335;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.214949;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.335536;
   return 0.;
}
float electronFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.114166;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.156868;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.117958;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0453808;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0352557;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0566433;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.01179;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0100706;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0137488;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0263784;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0220817;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0182064;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0541985;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0440799;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0301305;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.123434;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0976225;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0480083;
   return 0.;
}
float electronAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.400121;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.310303;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.434709;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.134924;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0702889;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.176986;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0678136;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0732791;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.128595;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0858109;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0980725;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.120241;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.206207;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.215078;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.204377;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.190836;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.221418;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.338085;
   return 0.;
}
float electronQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.473083;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.567948;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.392926;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.158189;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.221255;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.158681;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0711147;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0829365;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.142449;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0724598;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.12906;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.147341;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0876978;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.158748;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.181927;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.234259;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.159783;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.279726;
   return 0.;
}
float electronQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0386555;
   if (pt>=15 && pt<20 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0605446;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0831535;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0285307;
   if (pt>=20 && pt<25 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0397752;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0324968;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0123077;
   if (pt>=25 && pt<35 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0186246;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.023666;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0149396;
   if (pt>=35 && pt<50 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.030675;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0288891;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0338453;
   if (pt>=50 && pt<70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0478536;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0400865;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<0.8 ) return 0.0479921;
   if (pt>=70 && fabs(eta)>=0.8 && fabs(eta)<1.479 ) return 0.0422983;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.5 ) return 0.0376069;
   return 0.;
}
float muonFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.456206;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.503331;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.753583;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.134961;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.201255;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.240394;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0547376;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0825245;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.116965;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0320002;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.052804;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0977246;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0326371;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.079594;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.170545;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00868995;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.103198;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.136619;
   return 0.;
}
float muonFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0663701;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0774233;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.135457;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0231227;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0318872;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0620985;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00367513;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00574265;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0144333;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00529746;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00751073;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.02035;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0173518;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0193193;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0532072;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0517601;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0634385;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.124799;
   return 0.;
}
float muonAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.460885;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.506723;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.754469;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.13898;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.204947;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.24354;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0595521;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.086984;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.121878;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0461138;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0672001;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.115425;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0751219;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.111782;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.206497;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.118613;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.192278;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.220502;
   return 0.;
}
float muonQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.487579;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.546476;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.623434;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.180429;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.23425;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.277212;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0561123;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0862896;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0876005;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.045945;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0490797;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.062571;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0407374;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0527785;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0482802;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.042955;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0332809;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0815757;
   return 0.;
}
float muonQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0;
   if (pt>=10 && pt<15 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0;
   if (pt>=15 && pt<20 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0286323;
   if (pt>=15 && pt<20 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0377424;
   if (pt>=15 && pt<20 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0656204;
   if (pt>=20 && pt<25 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.0106157;
   if (pt>=20 && pt<25 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.0136068;
   if (pt>=20 && pt<25 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0278158;
   if (pt>=25 && pt<35 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00343789;
   if (pt>=25 && pt<35 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00629057;
   if (pt>=25 && pt<35 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0106443;
   if (pt>=35 && pt<50 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00378486;
   if (pt>=35 && pt<50 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00447896;
   if (pt>=35 && pt<50 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.011141;
   if (pt>=50 && pt<70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00441603;
   if (pt>=50 && pt<70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00874931;
   if (pt>=50 && pt<70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0108535;
   if (pt>=70 && fabs(eta)>=0 && fabs(eta)<1.2 ) return 0.00586428;
   if (pt>=70 && fabs(eta)>=1.2 && fabs(eta)<2.1 ) return 0.00641713;
   if (pt>=70 && fabs(eta)>=2.1 && fabs(eta)<2.4 ) return 0.0309623;
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

