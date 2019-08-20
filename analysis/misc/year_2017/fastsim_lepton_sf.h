float electronScaleFactorFastSim(float pt, float eta) {
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9551;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9667;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8783;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.7942;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.5977;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9296;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9297;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8320;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.8110;
  if (pt >= 20 && pt < 35 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.6604;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9345;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9314;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8175;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.8472;
  if (pt >= 35 && pt < 50 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.7485;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9374;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9356;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.7713;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.8546;
  if (pt >= 50 && pt < 100 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.7626;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9343;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9107;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.8821;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.8392;
  if (pt >= 100 && pt < 200 && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.8109;
  if (pt >= 200  && fabs(eta) >= 0.000 && fabs(eta) < 0.800) return 0.9659;
  if (pt >= 200  && fabs(eta) >= 0.800 && fabs(eta) < 1.444) return 0.9439;
  if (pt >= 200  && fabs(eta) >= 1.444 && fabs(eta) < 1.566) return 0.7670;
  if (pt >= 200  && fabs(eta) >= 1.566 && fabs(eta) < 2.000) return 0.7885;
  if (pt >= 200  && fabs(eta) >= 2.000 && fabs(eta) < 2.500) return 0.7022;
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
