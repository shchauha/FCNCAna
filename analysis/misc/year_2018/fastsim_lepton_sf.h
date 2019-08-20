float electronScaleFactorFastSim(float pt, float eta) {
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9333;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9145;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.7665;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.7098;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.6271;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9223;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9140;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8052;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.8242;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.7323;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9321;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9207;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8254;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.8691;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.7887;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9424;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9308;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8119;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.8817;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.8038;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9394;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9088;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.9107;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.9032;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.8629;
  if (pt >= 200  && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9659;
  if (pt >= 200  && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9406;
  if (pt >= 200  && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8443;
  if (pt >= 200  && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.7762;
  if (pt >= 200  && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.7272;
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
  return electronScaleFactorFastSim(pt, eta);
}
