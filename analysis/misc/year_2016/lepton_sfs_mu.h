
float muonScaleFactor(float pt, float eta) {
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.8964;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.8871;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9727;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9059;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9560;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9518;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9736;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9197;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9732;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9622;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9794;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9318;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9833;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9720;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9839;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9388;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9837;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9769;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9855;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9448;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9672;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9727;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9801;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9355;
  if (pt >= 60  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.9812;
  if (pt >= 60  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.9797;
  if (pt >= 60  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.9920;
  if (pt >= 60  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.9481;
  return 0.;
}

float muonScaleFactor_error(float pt, float eta) {
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0048;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0058;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0033;
  if (pt >= 10 && pt < 20 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0049;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0018;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0027;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0015;
  if (pt >= 20 && pt < 25 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0026;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0009;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0015;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0009;
  if (pt >= 25 && pt < 30 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0017;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0003;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0005;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0004;
  if (pt >= 30 && pt < 40 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0008;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0002;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0003;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0002;
  if (pt >= 40 && pt < 50 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0008;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0005;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0008;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0006;
  if (pt >= 50 && pt < 60 && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0019;
  if (pt >= 60  && fabs(eta) >= 0.000 && fabs(eta) < 0.900) return 0.0103;
  if (pt >= 60  && fabs(eta) >= 0.900 && fabs(eta) < 1.200) return 0.0089;
  if (pt >= 60  && fabs(eta) >= 1.200 && fabs(eta) < 2.100) return 0.0016;
  if (pt >= 60  && fabs(eta) >= 2.100 && fabs(eta) < 2.400) return 0.0441;
  return 1.;
}

float trackingSF(float eta) {
    if (eta >= -2.400 && eta < -2.100) return 0.991237;
    if (eta >= -2.100 && eta < -1.600) return 0.994853;
    if (eta >= -1.600 && eta < -1.200) return 0.996413;
    if (eta >= -1.200 && eta < -0.900) return 0.997157;
    if (eta >= -0.900 && eta < -0.600) return 0.997512;
    if (eta >= -0.600 && eta < -0.300) return 0.99756;
    if (eta >= -0.300 && eta < -0.200) return 0.996745;
    if (eta >= -0.200 && eta < 0.200) return 0.996996;
    if (eta >= 0.200 && eta < 0.300) return 0.99772;
    if (eta >= 0.300 && eta < 0.600) return 0.998604;
    if (eta >= 0.600 && eta < 0.900) return 0.998321;
    if (eta >= 0.900 && eta < 1.200) return 0.997682;
    if (eta >= 1.200 && eta < 1.600) return 0.995252;
    if (eta >= 1.600 && eta < 2.100) return 0.994919;
    if (eta >= 2.100 && eta < 2.400) return 0.987334;
    return 1.0;
}
