float flipRate(float pt, float eta) {
    float scale = 1.41;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.91502811928e-05*scale;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000114537784561*scale;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000737718492991*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 2.43389626268e-05*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000120545701142*scale;
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000911854481931*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 4.00203445913e-05*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000201819989855*scale;
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00131322856994*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 5.04595724607e-05*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000310283223174*scale;
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00203571326446*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 4.47847086164e-05*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000445696399346*scale;
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00237298940033*scale;
   if (pt>=200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 9.1849290275e-05*scale;
   if (pt>=200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000662764740024*scale;
   if (pt>=200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.00329508417923*scale;
   return 0.;
}
float flipRateError(float pt, float eta) {
    float scale = 1.41;
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 3.54576253168e-06*scale; // 12.2%
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 8.57622206089e-06*scale; // 7.5%
   if (pt>=15.0 && pt<40.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 2.61956706819e-05*scale; // 3.6%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 3.19666759313e-06*scale; // 13.1%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 8.97151384071e-06*scale; // 7.4%
   if (pt>=40.0 && pt<60.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 2.7954941261e-05*scale; // 3.1%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 9.8650423099e-06*scale; // 24.7%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 2.8172408285e-05*scale; // 14.0%
   if (pt>=60.0 && pt<80.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 8.83883929238e-05*scale; // 6.7%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 1.63169530122e-05*scale; // 32.3%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 5.20977373224e-05*scale; // 16.8%
   if (pt>=80.0 && pt<100.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000168299550426*scale; // 8.3%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 1.25207275411e-05*scale; // 28.0%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 5.22901640093e-05*scale; // 11.7%
   if (pt>=100.0 && pt<200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000158180857166*scale; // 6.7%
   if (pt>=200.0 && fabs(eta)>=0.000 && fabs(eta)<0.800 ) return 4.29836660284e-05*scale; // 46.8%
   if (pt>=200.0 && fabs(eta)>=0.800 && fabs(eta)<1.479 ) return 0.000166353866483*scale; // 25.1%
   if (pt>=200.0 && fabs(eta)>=1.479 && fabs(eta)<2.500 ) return 0.000542152139961*scale; // 16.5%
   return 0.;
}
