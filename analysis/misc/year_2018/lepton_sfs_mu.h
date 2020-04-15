float muonScaleFactor_Medium(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9916;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 1.0018;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 1.0031;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9889;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9951;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9962;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9935;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9733;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 1.0004;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9994;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9981;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9786;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9980;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9971;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9960;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9762;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9965;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9945;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9939;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9720;
  if (pt >= 60  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9989;
  if (pt >= 60  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9985;
  if (pt >= 60  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9957;
  if (pt >= 60  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9806;
  return 1.0;
}

float muonScaleFactorError_Medium(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0019;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0032;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0016;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0035;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0009;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0014;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0007;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0025;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0003;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0004;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0004;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0007;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0002;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0003;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0002;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0009;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0005;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0009;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0009;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0024;
  if (pt >= 60  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0010;
  if (pt >= 60  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0020;
  if (pt >= 60  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0015;
  if (pt >= 60  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0064;
  return 0.02;
}

float muonScaleFactor_RunABCD(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9824;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9784;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 1.0153;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0511;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9913;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9855;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 1.0110;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0271;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9948;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9906;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 1.0042;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0103;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9960;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9949;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 1.0010;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0041;
  if (pt >= 50  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9990;
  if (pt >= 50  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9970;
  if (pt >= 50  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 1.0010;
  if (pt >= 50  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 1.0030;
  return 1.0;
}

float muonScaleFactorError_RunABCD(float pt, float eta) {
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0224;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0157;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0255;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0106;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0019;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0027;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0068;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0060;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0015;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0018;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0015;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0027;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0014;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0014;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0020;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0014;
  if (pt >= 50  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0017;
  if (pt >= 50  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0025;
  if (pt >= 50  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0020;
  if (pt >= 50  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0027;
  return 0.02;
}

