float electronScaleFactorFastSimHighHT(float pt, float eta) {
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8448;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7726;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.5267;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8587;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7216;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.5905;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8637;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7067;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.6583;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8704;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7172;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.6979;
  if (pt >= 50  && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8662;
  if (pt >= 50  && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7414;
  if (pt >= 50  && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.7433;
  return 1.;
}

float electronScaleFactorFastSimLowHT(float pt, float eta) {
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8476;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7769;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.5291;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8574;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7218;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.5907;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8631;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7058;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.6569;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8708;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7196;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.6976;
  if (pt >= 50  && fabs(eta) >= 0.000 && fabs(eta) < 1.442) return 0.8724;
  if (pt >= 50  && fabs(eta) >= 1.442 && fabs(eta) < 1.556) return 0.7500;
  if (pt >= 50  && fabs(eta) >= 1.556 && fabs(eta) < 2.500) return 0.7451;
  return 1.;
}
float muonScaleFactorFastSim(float pt, float eta) {
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9575;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9191;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.8769;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9082;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9728;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9499;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9055;
  if (pt >= 20 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9513;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9826;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9693;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9318;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0154;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9876;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9805;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9455;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0466;
  if (pt >= 50  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9877;
  if (pt >= 50  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9821;
  if (pt >= 50  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9591;
  if (pt >= 50  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0935;
  return 1.;
}

float leptonScaleFactorFastSim(int pdgId, float pt, float eta, int ht) {
  if (abs(pdgId)==13) return muonScaleFactorFastSim(pt, eta);
  else if (abs(pdgId)==11){
    if (ht>300) return electronScaleFactorFastSimHighHT(pt, eta);
      else return electronScaleFactorFastSimLowHT(pt, eta);
  }
  return 0.;
}
