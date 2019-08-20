float flipRate(float pt, float eta) {
    float scale = 1.44;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.38453587652e-05*scale;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 7.99085449374e-05*scale;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000712067438702*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 1.91659933248e-05*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.0001166188839*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000779650545954*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.55933306063e-05*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000175292163676*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00130593692089*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 3.36222207333e-05*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000231082811867*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.0017362496487*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 5.82240178793e-05*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000378059688439*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00257346989023*scale;
   if (pt>=200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 0.000242357736611*scale;
   if (pt>=200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000816017591267*scale;
   if (pt>=200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00463215704213*scale;
   return 0.;
}
float flipRateError(float pt, float eta) {
    float scale = 1.44;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 3.21359412871e-06*scale; // 13.5%
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 7.33356108428e-06*scale; // 9.2%
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 2.51028034341e-05*scale; // 3.5%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.84478887391e-06*scale; // 14.8%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 8.66389018595e-06*scale; // 7.4%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 2.60352194351e-05*scale; // 3.3%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 7.78728063616e-06*scale; // 30.4%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 2.72187879158e-05*scale; // 15.5%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 9.34395924239e-05*scale; // 7.2%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 1.23509058561e-05*scale; // 36.7%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 4.88664968285e-05*scale; // 21.1%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000174139901881*scale; // 10.0%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 1.61381882453e-05*scale; // 27.7%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 5.79278577187e-05*scale; // 15.3%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000211062431866*scale; // 8.2%
   if (pt>=200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 0.000109698014086*scale; // 45.3%
   if (pt>=200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000269028749398*scale; // 33.0%
   if (pt>=200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00104557819978*scale; // 22.6%
   return 0.;
}
