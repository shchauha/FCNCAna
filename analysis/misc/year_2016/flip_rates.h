float flipRate(float pt, float eta) {
    float scale = 1.01;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 3.90282474065e-05*scale;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000196650768861*scale;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00131720111326*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.81181639943e-05*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000261378549702*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00175887162593*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.74709845052e-05*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000379026939692*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00304983346467*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 0.000161324430589*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000344730687906*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00335785460217*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 8.31256456405e-05*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.00076993412949*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00488458183303*scale;
   if (pt>=200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 7.56761225094e-05*scale;
   if (pt>=200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000315366076598*scale;
   if (pt>=200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00606167770218*scale;
   return 0.;
}
float flipRateError(float pt, float eta) {
    float scale = 1.01;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 5.83282386202e-06*scale; // 14.9%
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 1.59120785765e-05*scale; // 8.1%
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 4.9222388261e-05*scale; // 3.7%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 4.89371940015e-06*scale; // 17.4%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 1.8425858063e-05*scale; // 7.0%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 6.09206046183e-05*scale; // 3.5%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 1.12282455819e-05*scale; // 40.9%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 5.28678440254e-05*scale; // 13.9%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000199311523969*scale; // 6.5%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 4.0338362961e-05*scale; // 25.0%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 7.30127430122e-05*scale; // 21.2%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000315417639146*scale; // 9.4%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.22148082771e-05*scale; // 26.7%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 8.65314211887e-05*scale; // 11.2%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000305123823464*scale; // 6.2%
   if (pt>=200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 5.3575563497e-05*scale; // 70.8%
   if (pt>=200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000141773154432*scale; // 45.0%
   if (pt>=200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000947054283257*scale; // 15.6%
   return 0.;
}
