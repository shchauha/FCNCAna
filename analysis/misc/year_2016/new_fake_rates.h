float electronFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2301; // +-6.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2029; // +-6.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1861; // +-9.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0746; // +-10.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0942; // +-8.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0756; // +-11.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0647; // +-4.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0722; // +-3.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0707; // +-3.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0657; // +-7.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0739; // +-6.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0699; // +-5.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1117; // +-12.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1035; // +-11.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1135; // +-6.8%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2709; // +-13.7%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2297; // +-12.7%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2105; // +-8.1%;
   return 0.;
}
float electronFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0146;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0139;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0173;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0076;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0080;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0084;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0027;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0027;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0026;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0051;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0045;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0035;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0139;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0115;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0077;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0372;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0292;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0170;
   return 0.;
}

float electronAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2325; // +-6.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2046; // +-6.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1876; // +-9.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0775; // +-9.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0956; // +-8.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0771; // +-10.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0705; // +-3.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0755; // +-3.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0731; // +-3.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0884; // +-5.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0900; // +-4.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0796; // +-4.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1481; // +-8.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1300; // +-8.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1267; // +-6.0%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3380; // +-9.8%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2763; // +-9.8%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2312; // +-7.1%;
   return 0.;
}
float electronAlternativeFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0145;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0139;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0172;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0076;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0080;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0084;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0027;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0027;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0026;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0050;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0044;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0035;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0131;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0111;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0076;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0332;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0271;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0164;
   return 0.;
}

float electronQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2836; // +-3.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2852; // +-4.5%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2732; // +-8.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0883; // +-7.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0969; // +-9.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0777; // +-16.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0511; // +-11.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0643; // +-12.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0559; // +-16.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0906; // +-17.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0800; // +-14.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0804; // +-13.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1529; // +-15.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0917; // +-18.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0941; // +-14.8%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2224; // +-18.1%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1658; // +-11.2%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2602; // +-11.5%;
   return 0.;
}
float electronQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0112;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0127;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0244;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0069;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0095;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0129;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0059;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0080;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0091;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0160;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0114;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0108;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0243;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0165;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0140;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0402;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0186;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0298;
   return 0.;
}

float muonFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2916; // +-3.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3517; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3643; // +-7.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0899; // +-6.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0824; // +-8.5%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0847; // +-19.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0450; // +-10.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0514; // +-13.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0659; // +-25.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0432; // +-2.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0382; // +-3.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0413; // +-7.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0452; // +-4.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0367; // +-5.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0373; // +-14.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0455; // +-14.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0460; // +-15.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0194; // +-64.9%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0982; // +-20.1%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0838; // +-33.9%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0701; // +-110.1%;
   return 0.;
}
float muonFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0093;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0128;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0278;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0057;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0070;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0161;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0049;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0067;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0165;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0012;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0013;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0031;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0020;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0021;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0054;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0065;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0070;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0126;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0198;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0284;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0772;
   return 0.;
}

float muonAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2919; // +-3.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3519; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3645; // +-7.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0903; // +-6.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0827; // +-8.5%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0851; // +-18.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0456; // +-10.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0520; // +-12.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0665; // +-24.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0447; // +-2.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0394; // +-3.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0431; // +-7.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0512; // +-4.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0420; // +-5.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0452; // +-12.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0678; // +-9.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0646; // +-11.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0388; // +-33.2%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1458; // +-15.0%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1475; // +-19.7%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1575; // +-41.1%;
   return 0.;
}
float muonAlternativeFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0093;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0128;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0277;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0057;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0070;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0161;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0049;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0067;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0165;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0012;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0013;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0032;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0021;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0022;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0056;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0067;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0073;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0129;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0218;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0291;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0647;
   return 0.;
}

float muonQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3553; // +-1.1%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3758; // +-1.4%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3959; // +-2.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1004; // +-1.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1029; // +-1.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1165; // +-3.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0600; // +-2.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0580; // +-2.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0692; // +-4.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0513; // +-2.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0460; // +-3.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0528; // +-5.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0494; // +-4.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0414; // +-5.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0493; // +-9.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0483; // +-8.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0388; // +-10.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0474; // +-19.3%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0431; // +-17.6%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0487; // +-18.9%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0363; // +-43.9%;
   return 0.;
}
float muonQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0040;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0053;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0104;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0014;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0019;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0038;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0012;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0015;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0033;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0012;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0014;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0029;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0020;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0021;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0045;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0042;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0042;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0092;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0076;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0092;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0160;
   return 0.;
}

float electronFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3073; // +-5.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3143; // +-5.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2106; // +-7.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1633; // +-7.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1831; // +-6.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1246; // +-8.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1479; // +-4.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1542; // +-3.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1222; // +-3.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1717; // +-7.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1773; // +-5.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1452; // +-4.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2748; // +-10.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2312; // +-10.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2009; // +-6.5%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.5208; // +-9.5%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3851; // +-11.3%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3043; // +-7.4%;
   return 0.;
}
float electronFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0162;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0171;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0157;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0120;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0122;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0105;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0061;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0058;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0044;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0124;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0101;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0068;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0293;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0241;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0131;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0496;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0437;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0226;
   return 0.;
}

float electronAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3105; // +-5.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3165; // +-5.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2123; // +-7.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1678; // +-7.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1856; // +-6.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1265; // +-8.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1605; // +-3.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1612; // +-3.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1261; // +-3.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2240; // +-5.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2126; // +-4.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1637; // +-4.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3448; // +-7.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2843; // +-7.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2231; // +-5.6%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.6000; // +-6.4%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.4520; // +-8.1%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3310; // +-6.5%;
   return 0.;
}
float electronAlternativeFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0161;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0170;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0157;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0119;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0121;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0104;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0060;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0057;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0044;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0112;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0095;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0066;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0252;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0215;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0125;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0381;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0367;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0214;
   return 0.;
}

float electronQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3592; // +-3.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3730; // +-4.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2980; // +-9.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1555; // +-8.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1568; // +-10.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1180; // +-16.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1271; // +-11.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1335; // +-13.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0843; // +-17.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2265; // +-15.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1694; // +-13.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1359; // +-13.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3245; // +-14.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1490; // +-19.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1545; // +-14.6%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.4396; // +-15.2%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2775; // +-12.1%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3243; // +-11.0%;
   return 0.;
}
float electronQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0138;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0155;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0280;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0126;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0160;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0190;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0150;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0178;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0143;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0359;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0233;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0183;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0470;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0283;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0226;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0670;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0335;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0358;
   return 0.;
}

float muonFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3386; // +-2.1%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3714; // +-2.5%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3594; // +-5.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1056; // +-4.5%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1043; // +-5.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1003; // +-12.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0739; // +-7.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0744; // +-9.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0551; // +-20.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0613; // +-2.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0581; // +-3.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0515; // +-7.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0689; // +-4.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0604; // +-5.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0475; // +-15.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0719; // +-14.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0800; // +-15.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0638; // +-51.2%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1740; // +-20.0%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1093; // +-38.3%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0617; // +-138.2%;
   return 0.;
}
float muonFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0071;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0094;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0209;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0048;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0061;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0128;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0052;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0067;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0115;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0015;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0018;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0036;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0031;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0035;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0074;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0106;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0126;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0326;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0349;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0418;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0853;
   return 0.;
}

float muonAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3389; // +-2.1%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3716; // +-2.5%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3598; // +-5.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1061; // +-4.5%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1048; // +-5.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1008; // +-12.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0749; // +-7.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0753; // +-8.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0560; // +-20.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0637; // +-2.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0603; // +-3.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0539; // +-6.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0790; // +-3.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0699; // +-5.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0593; // +-12.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1097; // +-9.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1156; // +-10.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1135; // +-25.7%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2670; // +-11.7%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.2058; // +-18.2%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1585; // +-49.2%;
   return 0.;
}
float muonAlternativeFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0071;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0094;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0209;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0048;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0061;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0128;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0052;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0067;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0115;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0015;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0019;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0037;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0031;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0036;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0075;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0102;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0121;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0291;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0312;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0375;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0779;
   return 0.;
}

float muonQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.3760; // +-1.1%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3948; // +-1.4%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.4077; // +-2.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1218; // +-1.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1237; // +-1.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1334; // +-3.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0822; // +-2.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0777; // +-2.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0859; // +-4.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0772; // +-2.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0672; // +-3.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0721; // +-5.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0762; // +-3.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0643; // +-5.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0712; // +-9.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0798; // +-8.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0654; // +-10.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0720; // +-19.1%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0777; // +-17.3%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0811; // +-18.7%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0637; // +-43.2%;
   return 0.;
}
float muonQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0042;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0056;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0107;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0017;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0023;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0044;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0017;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0021;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0041;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0020;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0030;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0032;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0065;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0068;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0071;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0137;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0134;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0152;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0275;
   return 0.;
}

float fakeRate(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
    // SS uses iso for 2016; FT still uses nonisoHT trigs for ht>300 in 2016
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronFakeRate(pt,eta);
    else if (abs(id)==13) return muonFakeRate(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronFakeRate(pt,eta) : electronFakeRate_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonFakeRate(pt,eta) : muonFakeRate_IsoTrigs(pt,eta));
    else return 0.;
  }
}

float fakeRateError(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronFakeRateError(pt,eta);
    else if (abs(id)==13) return muonFakeRateError(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronFakeRateError(pt,eta) : electronFakeRateError_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonFakeRateError(pt,eta) : muonFakeRateError_IsoTrigs(pt,eta));
    else return 0.;
  }
}

float alternativeFakeRate(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronAlternativeFakeRate(pt,eta);
    else if (abs(id)==13) return muonAlternativeFakeRate(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronAlternativeFakeRate(pt,eta) : electronAlternativeFakeRate_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonAlternativeFakeRate(pt,eta) : muonAlternativeFakeRate_IsoTrigs(pt,eta));
    else return 0.;
  }
}

float qcdMCFakeRate(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
  if (ht>300. and ana != SSANA) {
    if (abs(id)==11) return electronQCDMCFakeRate(pt,eta);
    else if (abs(id)==13) return muonQCDMCFakeRate(pt,eta);
    else return 0.;
  } else {
    if (abs(id)==11) return (isLL ? electronQCDMCFakeRate(pt,eta) : electronQCDMCFakeRate_IsoTrigs(pt,eta));
    else if (abs(id)==13) return (isLL ? muonQCDMCFakeRate(pt,eta) : muonQCDMCFakeRate_IsoTrigs(pt,eta));
    else return 0.;
  }
}

float qcdMCFakeRateError(int id, float pt, float eta, float ht, ana_t ana, bool isLL) { 
    return 0.01;
}

