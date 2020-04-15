float electronFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1079; // +-19.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1445; // +-19.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1514; // +-16.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0814; // +-24.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0361; // +-32.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0960; // +-35.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0235; // +-21.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0307; // +-13.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0448; // +-8.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0235; // +-21.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0107; // +-78.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0449; // +-13.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0239; // +-114.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0199; // +-92.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0609; // +-19.7%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0239; // +-114.0%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0199; // +-92.8%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1061; // +-24.6%;
   return 0.;
}
float electronFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0206;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0277;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0251;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0195;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0117;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0340;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0051;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0041;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0051;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0084;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0058;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0273;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0185;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0120;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0273;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0185;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0261;
   return 0.;
}

float electronAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1126; // +-18.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1480; // +-18.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1539; // +-16.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0859; // +-22.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0398; // +-29.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0978; // +-34.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0334; // +-15.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0373; // +-11.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0486; // +-8.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0233; // +-33.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0420; // +-18.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0606; // +-9.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0939; // +-25.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0616; // +-28.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0832; // +-13.9%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1248; // +-52.6%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1232; // +-43.1%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1436; // +-16.9%;
   return 0.;
}
float electronAlternativeFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0205;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0275;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0250;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0194;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0116;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0338;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0051;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0041;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0079;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0079;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0058;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0237;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0172;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0116;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0656;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0531;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0243;
   return 0.;
}

float electronQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2035; // +-5.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2369; // +-6.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2406; // +-7.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0511; // +-8.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0711; // +-9.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0761; // +-10.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0373; // +-9.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0462; // +-10.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0489; // +-11.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0416; // +-15.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0472; // +-16.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0507; // +-15.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0446; // +-19.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0981; // +-19.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0780; // +-17.5%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1370; // +-18.7%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1299; // +-21.2%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1798; // +-15.4%;
   return 0.;
}
float electronQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0105;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0142;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0168;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0044;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0065;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0077;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0035;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0050;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0056;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0064;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0077;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0077;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0087;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0194;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0136;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0256;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0275;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0277;
   return 0.;
}

float muonFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2862; // +-3.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3393; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2842; // +-9.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0737; // +-6.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0628; // +-9.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0543; // +-19.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0401; // +-11.6%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0339; // +-14.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0215; // +-29.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0327; // +-5.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0303; // +-6.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0235; // +-13.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0323; // +-12.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0245; // +-16.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0353; // +-24.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0096; // +-137.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0245; // +-57.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0696; // +-52.3%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0096; // +-137.5%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1130; // +-36.9%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0696; // +-52.3%;
   return 0.;
}
float muonFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0091;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0122;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0266;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0050;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0059;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0104;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0047;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0050;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0064;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0020;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0032;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0040;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0087;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0132;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0141;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0364;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0132;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0417;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0364;
   return 0.;
}

float muonAlternativeFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2867; // +-3.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3396; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2844; // +-9.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0744; // +-6.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0633; // +-9.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0547; // +-18.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0412; // +-11.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0347; // +-14.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0223; // +-28.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0354; // +-5.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0325; // +-6.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0250; // +-12.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0431; // +-9.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0328; // +-11.9%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0441; // +-19.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0468; // +-26.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0545; // +-24.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1036; // +-33.1%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1027; // +-36.0%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1792; // +-20.9%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0297; // +-198.7%;
   return 0.;
}
float muonAlternativeFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0091;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0122;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0266;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0049;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0059;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0103;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0047;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0050;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0064;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0018;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0020;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0032;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0039;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0087;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0124;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0135;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0343;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0370;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0375;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0590;
   return 0.;
}

float muonQCDMCFakeRate(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2972; // +-1.8%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3197; // +-2.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2493; // +-6.1%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0809; // +-3.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0663; // +-4.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0576; // +-11.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0490; // +-5.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0368; // +-7.2%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0324; // +-14.9%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0366; // +-5.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0315; // +-7.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0315; // +-13.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0356; // +-8.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0313; // +-10.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0217; // +-25.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0495; // +-18.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0296; // +-28.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0438; // +-62.7%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0208; // +-42.1%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0244; // +-85.1%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0517; // +-66.7%;
   return 0.;
}
float muonQCDMCFakeRateError(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0052;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0073;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0152;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0029;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0032;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0068;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0025;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0027;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0048;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0021;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0024;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0043;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0031;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0033;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0056;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0092;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0083;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0275;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0088;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0208;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0345;
   return 0.;
}

float electronFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1959; // +-6.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2040; // +-7.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2122; // +-7.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0879; // +-9.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0807; // +-10.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1237; // +-7.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0586; // +-18.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0692; // +-13.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0951; // +-8.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0586; // +-18.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0347; // +-69.8%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1116; // +-14.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0663; // +-104.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0545; // +-92.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1202; // +-21.9%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0663; // +-104.2%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0545; // +-92.0%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1410; // +-28.4%;
   return 0.;
}
float electronFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0133;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0148;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0149;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0086;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0088;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0094;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0107;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0095;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0078;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0107;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0242;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0162;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0691;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0501;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0264;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0691;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0501;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0400;
   return 0.;
}

float electronAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2016; // +-6.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2093; // +-7.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2155; // +-6.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0960; // +-8.8%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0871; // +-10.0%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1272; // +-7.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0849; // +-12.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0839; // +-11.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1027; // +-7.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0723; // +-34.5%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1123; // +-18.4%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1458; // +-10.4%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2121; // +-23.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1508; // +-27.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1643; // +-14.6%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2489; // +-52.6%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2483; // +-36.5%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1940; // +-18.5%;
   return 0.;
}
float electronAlternativeFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0132;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0146;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0147;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0084;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0087;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0093;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0102;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0093;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0077;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0249;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0207;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0151;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0508;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0411;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0241;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1310;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0905;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0359;
   return 0.;
}

float electronQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000; // +-0.0%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000; // +-0.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2628; // +-5.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.3022; // +-6.0%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2998; // +-7.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1029; // +-8.7%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1355; // +-9.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1360; // +-10.8%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1063; // +-9.6%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1029; // +-11.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1056; // +-12.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1163; // +-15.0%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.1121; // +-16.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1114; // +-15.8%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.1017; // +-21.6%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2120; // +-18.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.1307; // +-18.3%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.2247; // +-19.2%;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.2284; // +-21.3%;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.2190; // +-14.4%;
   return 0.;
}
float electronQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0000;
   if (pt>=10 && pt<15 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0000;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0139;
   if (pt>=15 && pt<20 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0182;
   if (pt>=15 && pt<20 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0218;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0089;
   if (pt>=20 && pt<25 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0126;
   if (pt>=20 && pt<25 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0147;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0102;
   if (pt>=25 && pt<35 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0115;
   if (pt>=25 && pt<35 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0129;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0174;
   if (pt>=35 && pt<50 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0188;
   if (pt>=35 && pt<50 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0176;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0220;
   if (pt>=50 && pt<70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0385;
   if (pt>=50 && pt<70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0240;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<0.800) return 0.0432;
   if (pt>=70 && fabs(eta)>=0.800 && fabs(eta)<1.479) return 0.0487;
   if (pt>=70 && fabs(eta)>=1.479 && fabs(eta)<2.500) return 0.0315;
   return 0.;
}

float muonFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2890; // +-3.2%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3414; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2641; // +-9.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0748; // +-6.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0654; // +-9.4%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0495; // +-19.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0428; // +-11.4%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0426; // +-14.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0327; // +-31.4%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0361; // +-5.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0373; // +-6.7%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0331; // +-14.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0394; // +-12.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0335; // +-16.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0492; // +-25.7%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0120; // +-138.1%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0313; // +-59.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0977; // +-53.4%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0120; // +-138.1%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.1407; // +-39.7%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0977; // +-53.4%;
   return 0.;
}
float muonFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0091;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0124;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0259;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0050;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0061;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0096;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0049;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0063;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0103;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0020;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0025;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0047;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0048;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0055;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0126;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0166;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0186;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0522;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0166;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0559;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0522;
   return 0.;
}

float muonAlternativeFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2895; // +-3.1%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3417; // +-3.6%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.2643; // +-9.8%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0756; // +-6.7%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0661; // +-9.3%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0500; // +-19.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0440; // +-11.1%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0436; // +-14.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0336; // +-30.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0391; // +-5.1%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0400; // +-6.3%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0354; // +-13.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0526; // +-9.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0451; // +-12.1%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0619; // +-20.5%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0585; // +-26.3%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0707; // +-24.9%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1436; // +-35.3%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.1289; // +-35.8%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.2243; // +-23.0%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0867; // +-185.7%;
   return 0.;
}
float muonAlternativeFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0091;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0124;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0259;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0050;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0061;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0096;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0049;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0063;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0102;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0020;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0025;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0047;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0048;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0054;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0127;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0154;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0176;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0506;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0462;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0515;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.1609;
   return 0.;
}

float muonQCDMCFakeRate_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.2951; // +-1.9%;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.3243; // +-2.3%;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.3620; // +-4.2%;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0762; // +-3.9%;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0878; // +-4.6%;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0916; // +-8.9%;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0491; // +-5.3%;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0563; // +-6.5%;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0446; // +-13.2%;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0409; // +-6.0%;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0440; // +-7.5%;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0463; // +-13.7%;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0435; // +-8.2%;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0401; // +-10.3%;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0591; // +-21.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0439; // +-15.2%;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0434; // +-17.0%;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0450; // +-21.8%;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0457; // +-12.5%;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0358; // +-16.7%;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0698; // +-38.2%;
   return 0.;
}
float muonQCDMCFakeRateError_IsoTrigs(float pt, float eta) {
   if (pt>=10 && pt<15 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0056;
   if (pt>=10 && pt<15 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0073;
   if (pt>=10 && pt<15 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0154;
   if (pt>=15 && pt<20 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0030;
   if (pt>=15 && pt<20 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0041;
   if (pt>=15 && pt<20 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0082;
   if (pt>=20 && pt<25 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0026;
   if (pt>=20 && pt<25 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0037;
   if (pt>=20 && pt<25 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0059;
   if (pt>=25 && pt<35 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0025;
   if (pt>=25 && pt<35 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0033;
   if (pt>=25 && pt<35 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0063;
   if (pt>=35 && pt<50 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0035;
   if (pt>=35 && pt<50 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0041;
   if (pt>=35 && pt<50 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0124;
   if (pt>=50 && pt<70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0067;
   if (pt>=50 && pt<70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0074;
   if (pt>=50 && pt<70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0098;
   if (pt>=70 && fabs(eta)>=0.000 && fabs(eta)<1.200) return 0.0057;
   if (pt>=70 && fabs(eta)>=1.200 && fabs(eta)<2.100) return 0.0060;
   if (pt>=70 && fabs(eta)>=2.100 && fabs(eta)<2.400) return 0.0266;
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

