float muonScaleFactor_Medium(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9946;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 1.0046;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9955;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9760;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9943;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9963;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9935;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9746;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9985;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9995;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9982;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9791;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9967;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9975;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9959;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9764;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9938;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9967;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9946;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9711;
  if (pt >= 60  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9973;
  if (pt >= 60  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9957;
  if (pt >= 60  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9964;
  if (pt >= 60  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9810;
  return 1.0;
}

float muonScaleFactorError_Medium(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0032;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0043;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0072;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0198;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0011;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0076;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0016;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0209;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0003;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0005;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0057;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0264;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0002;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0004;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0002;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0010;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0006;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0013;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0008;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0030;
  if (pt >= 60  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0013;
  if (pt >= 60  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0043;
  if (pt >= 60  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0020;
  if (pt >= 60  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0076;
  return 0.02;
}

float muonScaleFactor_RunBCDEF(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9920;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9932;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9945;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9904;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9957;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9913;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9957;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9958;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9969;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9948;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9979;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9980;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9980;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9970;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9990;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9990;
  if (pt >= 50  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9990;
  if (pt >= 50  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9990;
  if (pt >= 50  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9990;
  if (pt >= 50  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0030;
  return 1.0;
}

float muonScaleFactorError_RunBCDEF(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0071;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0162;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0062;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0066;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0056;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0081;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0067;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0044;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0015;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0018;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0018;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0018;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0018;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0018;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0017;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0017;
  if (pt >= 50  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0014;
  if (pt >= 50  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0014;
  if (pt >= 50  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0020;
  if (pt >= 50  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0020;
  return 0.02;
}

