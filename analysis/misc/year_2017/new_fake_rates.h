float electronFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1254; // +-21.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2051; // +-16.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1255; // +-23.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0636; // +-29.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0236; // +-37.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0774; // +-29.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0222; // +-17.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0322; // +-13.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0376; // +-11.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0113; // +-85.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0098; // +-82.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0302; // +-23.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0633; // +-37.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0588; // +-40.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0795; // +-19.4%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1141; // +-89.2%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1512; // +-39.5%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1510; // +-23.4%;
   return 0.;
}
float electronFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0273;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0330;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0296;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0185;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0087;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0229;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0038;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0043;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0044;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0097;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0081;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0071;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0240;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0235;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0154;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1018;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0597;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0353;
   return 0.;
}

float electronAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1291; // +-21.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2079; // +-15.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1286; // +-23.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0671; // +-27.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0263; // +-33.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0800; // +-28.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0304; // +-12.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0374; // +-11.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0421; // +-10.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0448; // +-20.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0353; // +-22.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0472; // +-14.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1071; // +-21.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0951; // +-23.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1012; // +-14.9%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2733; // +-25.9%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2255; // +-22.6%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1889; // +-17.9%;
   return 0.;
}
float electronAlternativeFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0272;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0329;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0295;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0184;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0087;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0227;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0039;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0043;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0044;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0092;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0078;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0070;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0227;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0222;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0151;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0709;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0510;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0338;
   return 0.;
}

float electronQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2133; // +-5.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2425; // +-5.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2442; // +-7.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0548; // +-9.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0803; // +-10.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0752; // +-11.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0309; // +-9.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0414; // +-11.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0563; // +-12.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0509; // +-18.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0457; // +-19.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0585; // +-17.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0458; // +-31.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0933; // +-26.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0824; // +-23.0%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1757; // +-18.4%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1353; // +-22.0%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2190; // +-14.6%;
   return 0.;
}
float electronQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0108;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0140;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0176;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0050;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0085;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0084;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0031;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0049;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0067;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0093;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0089;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0099;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0144;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0250;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0190;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0324;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0298;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0319;
   return 0.;
}

float muonFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2969; // +-4.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3285; // +-5.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2973; // +-14.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0679; // +-17.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0780; // +-11.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0664; // +-25.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0277; // +-18.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0284; // +-22.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0215; // +-51.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0309; // +-5.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0328; // +-6.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0365; // +-11.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0261; // +-13.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0270; // +-14.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0266; // +-30.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0399; // +-28.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0401; // +-28.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0610; // +-45.6%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0589; // +-61.3%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0538; // +-66.1%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0610; // +-45.6%;
   return 0.;
}
float muonFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0136;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0175;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0433;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0121;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0088;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0166;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0051;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0063;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0109;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0016;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0022;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0042;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0035;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0080;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0113;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0114;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0278;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0361;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0355;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0278;
   return 0.;
}

float muonAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2973; // +-4.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3288; // +-5.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2976; // +-14.5%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0683; // +-17.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0785; // +-11.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0668; // +-24.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0286; // +-18.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0291; // +-21.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0221; // +-49.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0330; // +-5.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0343; // +-6.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0381; // +-11.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0340; // +-10.6%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0344; // +-11.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0344; // +-23.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0641; // +-17.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0583; // +-19.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0811; // +-34.3%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1251; // +-28.6%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1029; // +-33.5%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0706; // +-106.6%;
   return 0.;
}
float muonAlternativeFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0135;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0175;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0433;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0121;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0088;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0166;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0051;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0063;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0109;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0017;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0023;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0042;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0036;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0082;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0113;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0116;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0278;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0358;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0344;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0753;
   return 0.;
}

float muonQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2950; // +-1.9%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3308; // +-2.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3517; // +-4.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0776; // +-4.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0893; // +-4.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0836; // +-9.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0505; // +-5.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0512; // +-6.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0398; // +-13.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0399; // +-5.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0370; // +-7.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0374; // +-13.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0421; // +-7.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0338; // +-9.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0456; // +-17.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0405; // +-11.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0402; // +-15.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0461; // +-21.0%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0473; // +-13.2%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0307; // +-19.1%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0776; // +-34.5%;
   return 0.;
}
float muonQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0057;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0076;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0153;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0031;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0043;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0077;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0028;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0034;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0054;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0023;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0029;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0049;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0032;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0031;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0080;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0047;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0063;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0097;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0063;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0059;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0268;
   return 0.;
}

float electronFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2409; // +-7.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2371; // +-8.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1796; // +-12.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1057; // +-11.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0906; // +-12.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1169; // +-12.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0618; // +-18.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0776; // +-13.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0867; // +-11.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0381; // +-79.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0327; // +-74.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0905; // +-20.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2528; // +-27.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1498; // +-36.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1930; // +-16.3%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2840; // +-54.3%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2858; // +-31.4%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2607; // +-19.2%;
   return 0.;
}
float electronFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0170;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0211;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0232;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0118;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0115;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0144;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0113;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0101;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0096;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0301;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0242;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0185;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0683;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0541;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0314;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1543;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0898;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0500;
   return 0.;
}

float electronAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2454; // +-6.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2414; // +-8.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1852; // +-12.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1122; // +-10.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0963; // +-11.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1216; // +-11.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0850; // +-12.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0899; // +-11.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0964; // +-9.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1291; // +-19.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1012; // +-20.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1303; // +-13.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3576; // +-14.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2245; // +-20.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2321; // +-12.5%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.4794; // +-18.4%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3821; // +-18.4%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3095; // +-14.7%;
   return 0.;
}
float electronAlternativeFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0168;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0209;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0229;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0116;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0114;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0142;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0108;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0099;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0094;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0249;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0211;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0170;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0530;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0463;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0289;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0884;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0704;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0455;
   return 0.;
}

float electronQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2823; // +-5.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3082; // +-5.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.3024; // +-7.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1045; // +-9.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1559; // +-10.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1180; // +-12.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0871; // +-10.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0966; // +-12.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1205; // +-12.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1197; // +-17.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1024; // +-20.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1246; // +-17.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0962; // +-34.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1901; // +-25.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1543; // +-23.2%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.3263; // +-17.2%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2156; // +-24.1%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2857; // +-14.5%;
   return 0.;
}
float electronQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0141;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0182;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0232;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0095;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0164;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0144;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0091;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0121;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0151;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0214;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0208;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0217;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0332;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0486;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0358;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0562;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0520;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0414;
   return 0.;
}

float muonFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2908; // +-4.4%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3353; // +-4.7%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2772; // +-12.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0798; // +-10.5%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0755; // +-11.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0942; // +-20.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0364; // +-15.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0406; // +-18.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0451; // +-30.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0319; // +-5.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0357; // +-6.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0365; // +-11.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0312; // +-13.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0291; // +-15.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0407; // +-30.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0365; // +-34.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0485; // +-31.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0523; // +-70.2%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0525; // +-78.6%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0508; // +-88.4%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0523; // +-70.2%;
   return 0.;
}
float muonFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0128;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0159;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0359;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0083;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0083;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0190;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0056;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0075;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0136;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0023;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0043;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0043;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0046;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0122;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0124;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0152;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0367;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0413;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0449;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0367;
   return 0.;
}

float muonAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2912; // +-4.4%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3356; // +-4.7%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2776; // +-12.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0803; // +-10.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0760; // +-10.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0947; // +-20.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0374; // +-14.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0414; // +-18.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0459; // +-29.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0342; // +-5.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0375; // +-6.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0383; // +-11.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0412; // +-10.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0376; // +-12.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0520; // +-23.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0638; // +-19.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0729; // +-20.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0825; // +-43.2%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1321; // +-28.4%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1131; // +-37.6%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1074; // +-154.1%;
   return 0.;
}
float muonAlternativeFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0127;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0159;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0359;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0083;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0083;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0190;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0056;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0075;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0136;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0023;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0044;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0043;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0046;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0123;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0122;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0150;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0357;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0375;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0425;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1655;
   return 0.;
}

float muonQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2951; // +-1.9%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3313; // +-2.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3530; // +-4.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0778; // +-4.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0902; // +-4.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0848; // +-9.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0509; // +-5.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0521; // +-6.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0405; // +-13.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0404; // +-5.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0381; // +-7.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0386; // +-13.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0428; // +-7.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0354; // +-9.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0478; // +-17.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0419; // +-11.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0427; // +-15.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0485; // +-21.0%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0491; // +-13.2%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0344; // +-18.9%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0834; // +-34.3%;
   return 0.;
}
float muonQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0057;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0077;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0154;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0031;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0043;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0078;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0028;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0035;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0055;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0023;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0029;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0051;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0033;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0033;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0083;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0048;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0066;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0102;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0065;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0065;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0286;
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

