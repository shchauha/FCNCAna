
float electronScaleFactorHighHT(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 0.8169;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.8139;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.442) return 0.9570;
  if (pt >= 10 && pt < 20 && eta >= -1.442 && eta < -0.800) return 0.9437;
  if (pt >= 10 && pt < 20 && eta >= -0.800 && eta < -0.000) return 0.9245;
  if (pt >= 10 && pt < 20 && eta >= -0.000 && eta < 0.800) return 0.9245;
  if (pt >= 10 && pt < 20 && eta >= 0.800 && eta < 1.442) return 0.9437;
  if (pt >= 10 && pt < 20 && eta >= 1.442 && eta < 1.566) return 0.9570;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.8139;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 0.8169;
  if (pt >= 20 && pt < 30 && eta >= -2.500 && eta < -2.000) return 0.8477;
  if (pt >= 20 && pt < 30 && eta >= -2.000 && eta < -1.566) return 0.8890;
  if (pt >= 20 && pt < 30 && eta >= -1.566 && eta < -1.442) return 0.9346;
  if (pt >= 20 && pt < 30 && eta >= -1.442 && eta < -0.800) return 0.9389;
  if (pt >= 20 && pt < 30 && eta >= -0.800 && eta < -0.000) return 0.9419;
  if (pt >= 20 && pt < 30 && eta >= -0.000 && eta < 0.800) return 0.9419;
  if (pt >= 20 && pt < 30 && eta >= 0.800 && eta < 1.442) return 0.9389;
  if (pt >= 20 && pt < 30 && eta >= 1.442 && eta < 1.566) return 0.9346;
  if (pt >= 20 && pt < 30 && eta >= 1.566 && eta < 2.000) return 0.8890;
  if (pt >= 20 && pt < 30 && eta >= 2.000 && eta < 2.500) return 0.8477;
  if (pt >= 30 && pt < 40 && eta >= -2.500 && eta < -2.000) return 0.9042;
  if (pt >= 30 && pt < 40 && eta >= -2.000 && eta < -1.566) return 0.9222;
  if (pt >= 30 && pt < 40 && eta >= -1.566 && eta < -1.442) return 0.9696;
  if (pt >= 30 && pt < 40 && eta >= -1.442 && eta < -0.800) return 0.9547;
  if (pt >= 30 && pt < 40 && eta >= -0.800 && eta < -0.000) return 0.9533;
  if (pt >= 30 && pt < 40 && eta >= -0.000 && eta < 0.800) return 0.9533;
  if (pt >= 30 && pt < 40 && eta >= 0.800 && eta < 1.442) return 0.9547;
  if (pt >= 30 && pt < 40 && eta >= 1.442 && eta < 1.566) return 0.9696;
  if (pt >= 30 && pt < 40 && eta >= 1.566 && eta < 2.000) return 0.9222;
  if (pt >= 30 && pt < 40 && eta >= 2.000 && eta < 2.500) return 0.9042;
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.000) return 0.9388;
  if (pt >= 40 && pt < 50 && eta >= -2.000 && eta < -1.566) return 0.9430;
  if (pt >= 40 && pt < 50 && eta >= -1.566 && eta < -1.442) return 0.9620;
  if (pt >= 40 && pt < 50 && eta >= -1.442 && eta < -0.800) return 0.9565;
  if (pt >= 40 && pt < 50 && eta >= -0.800 && eta < -0.000) return 0.9516;
  if (pt >= 40 && pt < 50 && eta >= -0.000 && eta < 0.800) return 0.9516;
  if (pt >= 40 && pt < 50 && eta >= 0.800 && eta < 1.442) return 0.9565;
  if (pt >= 40 && pt < 50 && eta >= 1.442 && eta < 1.566) return 0.9620;
  if (pt >= 40 && pt < 50 && eta >= 1.566 && eta < 2.000) return 0.9430;
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.500) return 0.9388;
  if (pt >= 50 && pt < 100 && eta >= -2.500 && eta < -2.000) return 0.9486;
  if (pt >= 50 && pt < 100 && eta >= -2.000 && eta < -1.566) return 0.9511;
  if (pt >= 50 && pt < 100 && eta >= -1.566 && eta < -1.442) return 0.9320;
  if (pt >= 50 && pt < 100 && eta >= -1.442 && eta < -0.800) return 0.9535;
  if (pt >= 50 && pt < 100 && eta >= -0.800 && eta < -0.000) return 0.9490;
  if (pt >= 50 && pt < 100 && eta >= -0.000 && eta < 0.800) return 0.9490;
  if (pt >= 50 && pt < 100 && eta >= 0.800 && eta < 1.442) return 0.9535;
  if (pt >= 50 && pt < 100 && eta >= 1.442 && eta < 1.566) return 0.9320;
  if (pt >= 50 && pt < 100 && eta >= 1.566 && eta < 2.000) return 0.9511;
  if (pt >= 50 && pt < 100 && eta >= 2.000 && eta < 2.500) return 0.9486;
  if (pt >= 100  && eta >= -2.500 && eta < -2.000) return 1.0145;
  if (pt >= 100  && eta >= -2.000 && eta < -1.566) return 0.9500;
  if (pt >= 100  && eta >= -1.566 && eta < -1.442) return 0.8957;
  if (pt >= 100  && eta >= -1.442 && eta < -0.800) return 0.9616;
  if (pt >= 100  && eta >= -0.800 && eta < -0.000) return 0.9573;
  if (pt >= 100  && eta >= -0.000 && eta < 0.800) return 0.9573;
  if (pt >= 100  && eta >= 0.800 && eta < 1.442) return 0.9616;
  if (pt >= 100  && eta >= 1.442 && eta < 1.566) return 0.8957;
  if (pt >= 100  && eta >= 1.566 && eta < 2.000) return 0.9500;
  if (pt >= 100  && eta >= 2.000 && eta < 2.500) return 1.0145;
  return 0.;
}

float electronScaleFactorLowHT(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 0.8165;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.8107;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.442) return 0.9568;
  if (pt >= 10 && pt < 20 && eta >= -1.442 && eta < -0.800) return 0.9437;
  if (pt >= 10 && pt < 20 && eta >= -0.800 && eta < -0.000) return 0.9212;
  if (pt >= 10 && pt < 20 && eta >= -0.000 && eta < 0.800) return 0.9212;
  if (pt >= 10 && pt < 20 && eta >= 0.800 && eta < 1.442) return 0.9437;
  if (pt >= 10 && pt < 20 && eta >= 1.442 && eta < 1.566) return 0.9568;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.8107;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 0.8165;
  if (pt >= 20 && pt < 30 && eta >= -2.500 && eta < -2.000) return 0.8455;
  if (pt >= 20 && pt < 30 && eta >= -2.000 && eta < -1.566) return 0.8897;
  if (pt >= 20 && pt < 30 && eta >= -1.566 && eta < -1.442) return 0.9357;
  if (pt >= 20 && pt < 30 && eta >= -1.442 && eta < -0.800) return 0.9385;
  if (pt >= 20 && pt < 30 && eta >= -0.800 && eta < -0.000) return 0.9401;
  if (pt >= 20 && pt < 30 && eta >= -0.000 && eta < 0.800) return 0.9401;
  if (pt >= 20 && pt < 30 && eta >= 0.800 && eta < 1.442) return 0.9385;
  if (pt >= 20 && pt < 30 && eta >= 1.442 && eta < 1.566) return 0.9357;
  if (pt >= 20 && pt < 30 && eta >= 1.566 && eta < 2.000) return 0.8897;
  if (pt >= 20 && pt < 30 && eta >= 2.000 && eta < 2.500) return 0.8455;
  if (pt >= 30 && pt < 40 && eta >= -2.500 && eta < -2.000) return 0.9044;
  if (pt >= 30 && pt < 40 && eta >= -2.000 && eta < -1.566) return 0.9224;
  if (pt >= 30 && pt < 40 && eta >= -1.566 && eta < -1.442) return 0.9701;
  if (pt >= 30 && pt < 40 && eta >= -1.442 && eta < -0.800) return 0.9543;
  if (pt >= 30 && pt < 40 && eta >= -0.800 && eta < -0.000) return 0.9522;
  if (pt >= 30 && pt < 40 && eta >= -0.000 && eta < 0.800) return 0.9522;
  if (pt >= 30 && pt < 40 && eta >= 0.800 && eta < 1.442) return 0.9543;
  if (pt >= 30 && pt < 40 && eta >= 1.442 && eta < 1.566) return 0.9701;
  if (pt >= 30 && pt < 40 && eta >= 1.566 && eta < 2.000) return 0.9224;
  if (pt >= 30 && pt < 40 && eta >= 2.000 && eta < 2.500) return 0.9044;
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.000) return 0.9389;
  if (pt >= 40 && pt < 50 && eta >= -2.000 && eta < -1.566) return 0.9431;
  if (pt >= 40 && pt < 50 && eta >= -1.566 && eta < -1.442) return 0.9616;
  if (pt >= 40 && pt < 50 && eta >= -1.442 && eta < -0.800) return 0.9562;
  if (pt >= 40 && pt < 50 && eta >= -0.800 && eta < -0.000) return 0.9506;
  if (pt >= 40 && pt < 50 && eta >= -0.000 && eta < 0.800) return 0.9506;
  if (pt >= 40 && pt < 50 && eta >= 0.800 && eta < 1.442) return 0.9562;
  if (pt >= 40 && pt < 50 && eta >= 1.442 && eta < 1.566) return 0.9616;
  if (pt >= 40 && pt < 50 && eta >= 1.566 && eta < 2.000) return 0.9431;
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.500) return 0.9389;
  if (pt >= 50 && pt < 100 && eta >= -2.500 && eta < -2.000) return 0.9488;
  if (pt >= 50 && pt < 100 && eta >= -2.000 && eta < -1.566) return 0.9516;
  if (pt >= 50 && pt < 100 && eta >= -1.566 && eta < -1.442) return 0.9324;
  if (pt >= 50 && pt < 100 && eta >= -1.442 && eta < -0.800) return 0.9532;
  if (pt >= 50 && pt < 100 && eta >= -0.800 && eta < -0.000) return 0.9484;
  if (pt >= 50 && pt < 100 && eta >= -0.000 && eta < 0.800) return 0.9484;
  if (pt >= 50 && pt < 100 && eta >= 0.800 && eta < 1.442) return 0.9532;
  if (pt >= 50 && pt < 100 && eta >= 1.442 && eta < 1.566) return 0.9324;
  if (pt >= 50 && pt < 100 && eta >= 1.566 && eta < 2.000) return 0.9516;
  if (pt >= 50 && pt < 100 && eta >= 2.000 && eta < 2.500) return 0.9488;
  if (pt >= 100  && eta >= -2.500 && eta < -2.000) return 1.0091;
  if (pt >= 100  && eta >= -2.000 && eta < -1.566) return 0.9358;
  if (pt >= 100  && eta >= -1.566 && eta < -1.442) return 0.8804;
  if (pt >= 100  && eta >= -1.442 && eta < -0.800) return 0.9390;
  if (pt >= 100  && eta >= -0.800 && eta < -0.000) return 0.9458;
  if (pt >= 100  && eta >= -0.000 && eta < 0.800) return 0.9458;
  if (pt >= 100  && eta >= 0.800 && eta < 1.442) return 0.9390;
  if (pt >= 100  && eta >= 1.442 && eta < 1.566) return 0.8804;
  if (pt >= 100  && eta >= 1.566 && eta < 2.000) return 0.9358;
  if (pt >= 100  && eta >= 2.000 && eta < 2.500) return 1.0091;
  return 0.;
}

float electronGSF(float pt, float eta) {
  if (pt >= 25  && eta >= -2.500 && eta < -2.450) return 1.3176;
  if (pt >= 25  && eta >= -2.450 && eta < -2.400) return 1.1138;
  if (pt >= 25  && eta >= -2.400 && eta < -2.300) return 1.0246;
  if (pt >= 25  && eta >= -2.300 && eta < -2.200) return 1.0136;
  if (pt >= 25  && eta >= -2.200 && eta < -2.000) return 1.0073;
  if (pt >= 25  && eta >= -2.000 && eta < -1.800) return 0.9948;
  if (pt >= 25  && eta >= -1.800 && eta < -1.630) return 0.9948;
  if (pt >= 25  && eta >= -1.630 && eta < -1.566) return 0.9916;
  if (pt >= 25  && eta >= -1.566 && eta < -1.444) return 0.9631;
  if (pt >= 25  && eta >= -1.444 && eta < -1.200) return 0.9897;
  if (pt >= 25  && eta >= -1.200 && eta < -1.000) return 0.9857;
  if (pt >= 25  && eta >= -1.000 && eta < -0.600) return 0.9816;
  if (pt >= 25  && eta >= -0.600 && eta < -0.400) return 0.9847;
  if (pt >= 25  && eta >= -0.400 && eta < -0.200) return 0.9816;
  if (pt >= 25  && eta >= -0.200 && eta < 0.000) return 0.9804;
  if (pt >= 25  && eta >= 0.000 && eta < 0.200) return 0.9846;
  if (pt >= 25  && eta >= 0.200 && eta < 0.400) return 0.9888;
  if (pt >= 25  && eta >= 0.400 && eta < 0.600) return 0.9877;
  if (pt >= 25  && eta >= 0.600 && eta < 1.000) return 0.9877;
  if (pt >= 25  && eta >= 1.000 && eta < 1.200) return 0.9877;
  if (pt >= 25  && eta >= 1.200 && eta < 1.444) return 0.9877;
  if (pt >= 25  && eta >= 1.444 && eta < 1.566) return 0.9676;
  if (pt >= 25  && eta >= 1.566 && eta < 1.630) return 0.9896;
  if (pt >= 25  && eta >= 1.630 && eta < 1.800) return 0.9928;
  if (pt >= 25  && eta >= 1.800 && eta < 2.000) return 0.9918;
  if (pt >= 25  && eta >= 2.000 && eta < 2.200) return 0.9979;
  if (pt >= 25  && eta >= 2.200 && eta < 2.300) return 1.0010;
  if (pt >= 25  && eta >= 2.300 && eta < 2.400) return 0.9895;
  if (pt >= 25  && eta >= 2.400 && eta < 2.450) return 0.9705;
  if (pt >= 25  && eta >= 2.450 && eta < 2.500) return 0.9067;
  return 1.0;
}


float electronScaleFactorHighHT_error(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 0.0316;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.0284;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.442) return 0.0639;
  if (pt >= 10 && pt < 20 && eta >= -1.442 && eta < -0.800) return 0.0119;
  if (pt >= 10 && pt < 20 && eta >= -0.800 && eta < -0.000) return 0.0355;
  if (pt >= 10 && pt < 20 && eta >= -0.000 && eta < 0.800) return 0.0355;
  if (pt >= 10 && pt < 20 && eta >= 0.800 && eta < 1.442) return 0.0119;
  if (pt >= 10 && pt < 20 && eta >= 1.442 && eta < 1.566) return 0.0639;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.0284;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 0.0316;
  if (pt >= 20 && pt < 30 && eta >= -2.500 && eta < -2.000) return 0.0177;
  if (pt >= 20 && pt < 30 && eta >= -2.000 && eta < -1.566) return 0.0173;
  if (pt >= 20 && pt < 30 && eta >= -1.566 && eta < -1.442) return 0.1184;
  if (pt >= 20 && pt < 30 && eta >= -1.442 && eta < -0.800) return 0.0258;
  if (pt >= 20 && pt < 30 && eta >= -0.800 && eta < -0.000) return 0.0169;
  if (pt >= 20 && pt < 30 && eta >= -0.000 && eta < 0.800) return 0.0169;
  if (pt >= 20 && pt < 30 && eta >= 0.800 && eta < 1.442) return 0.0258;
  if (pt >= 20 && pt < 30 && eta >= 1.442 && eta < 1.566) return 0.1184;
  if (pt >= 20 && pt < 30 && eta >= 1.566 && eta < 2.000) return 0.0173;
  if (pt >= 20 && pt < 30 && eta >= 2.000 && eta < 2.500) return 0.0177;
  if (pt >= 30 && pt < 40 && eta >= -2.500 && eta < -2.000) return 0.0090;
  if (pt >= 30 && pt < 40 && eta >= -2.000 && eta < -1.566) return 0.0068;
  if (pt >= 30 && pt < 40 && eta >= -1.566 && eta < -1.442) return 0.0089;
  if (pt >= 30 && pt < 40 && eta >= -1.442 && eta < -0.800) return 0.0084;
  if (pt >= 30 && pt < 40 && eta >= -0.800 && eta < -0.000) return 0.0070;
  if (pt >= 30 && pt < 40 && eta >= -0.000 && eta < 0.800) return 0.0070;
  if (pt >= 30 && pt < 40 && eta >= 0.800 && eta < 1.442) return 0.0084;
  if (pt >= 30 && pt < 40 && eta >= 1.442 && eta < 1.566) return 0.0089;
  if (pt >= 30 && pt < 40 && eta >= 1.566 && eta < 2.000) return 0.0068;
  if (pt >= 30 && pt < 40 && eta >= 2.000 && eta < 2.500) return 0.0090;
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.000) return 0.0146;
  if (pt >= 40 && pt < 50 && eta >= -2.000 && eta < -1.566) return 0.0066;
  if (pt >= 40 && pt < 50 && eta >= -1.566 && eta < -1.442) return 0.0168;
  if (pt >= 40 && pt < 50 && eta >= -1.442 && eta < -0.800) return 0.0040;
  if (pt >= 40 && pt < 50 && eta >= -0.800 && eta < -0.000) return 0.0036;
  if (pt >= 40 && pt < 50 && eta >= -0.000 && eta < 0.800) return 0.0036;
  if (pt >= 40 && pt < 50 && eta >= 0.800 && eta < 1.442) return 0.0040;
  if (pt >= 40 && pt < 50 && eta >= 1.442 && eta < 1.566) return 0.0168;
  if (pt >= 40 && pt < 50 && eta >= 1.566 && eta < 2.000) return 0.0066;
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.500) return 0.0146;
  if (pt >= 50 && pt < 100 && eta >= -2.500 && eta < -2.000) return 0.0197;
  if (pt >= 50 && pt < 100 && eta >= -2.000 && eta < -1.566) return 0.0110;
  if (pt >= 50 && pt < 100 && eta >= -1.566 && eta < -1.442) return 0.0171;
  if (pt >= 50 && pt < 100 && eta >= -1.442 && eta < -0.800) return 0.0056;
  if (pt >= 50 && pt < 100 && eta >= -0.800 && eta < -0.000) return 0.0084;
  if (pt >= 50 && pt < 100 && eta >= -0.000 && eta < 0.800) return 0.0084;
  if (pt >= 50 && pt < 100 && eta >= 0.800 && eta < 1.442) return 0.0056;
  if (pt >= 50 && pt < 100 && eta >= 1.442 && eta < 1.566) return 0.0171;
  if (pt >= 50 && pt < 100 && eta >= 1.566 && eta < 2.000) return 0.0110;
  if (pt >= 50 && pt < 100 && eta >= 2.000 && eta < 2.500) return 0.0197;
  if (pt >= 100  && eta >= -2.500 && eta < -2.000) return 0.0605;
  if (pt >= 100  && eta >= -2.000 && eta < -1.566) return 0.0145;
  if (pt >= 100  && eta >= -1.566 && eta < -1.442) return 0.0655;
  if (pt >= 100  && eta >= -1.442 && eta < -0.800) return 0.0081;
  if (pt >= 100  && eta >= -0.800 && eta < -0.000) return 0.0060;
  if (pt >= 100  && eta >= -0.000 && eta < 0.800) return 0.0060;
  if (pt >= 100  && eta >= 0.800 && eta < 1.442) return 0.0081;
  if (pt >= 100  && eta >= 1.442 && eta < 1.566) return 0.0655;
  if (pt >= 100  && eta >= 1.566 && eta < 2.000) return 0.0145;
  if (pt >= 100  && eta >= 2.000 && eta < 2.500) return 0.0605;
  return 1.;
}

float electronScaleFactorLowHT_error(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 0.0325;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.0288;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.442) return 0.0645;
  if (pt >= 10 && pt < 20 && eta >= -1.442 && eta < -0.800) return 0.0108;
  if (pt >= 10 && pt < 20 && eta >= -0.800 && eta < -0.000) return 0.0355;
  if (pt >= 10 && pt < 20 && eta >= -0.000 && eta < 0.800) return 0.0355;
  if (pt >= 10 && pt < 20 && eta >= 0.800 && eta < 1.442) return 0.0108;
  if (pt >= 10 && pt < 20 && eta >= 1.442 && eta < 1.566) return 0.0645;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.0288;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 0.0325;
  if (pt >= 20 && pt < 30 && eta >= -2.500 && eta < -2.000) return 0.0172;
  if (pt >= 20 && pt < 30 && eta >= -2.000 && eta < -1.566) return 0.0170;
  if (pt >= 20 && pt < 30 && eta >= -1.566 && eta < -1.442) return 0.1187;
  if (pt >= 20 && pt < 30 && eta >= -1.442 && eta < -0.800) return 0.0259;
  if (pt >= 20 && pt < 30 && eta >= -0.800 && eta < -0.000) return 0.0176;
  if (pt >= 20 && pt < 30 && eta >= -0.000 && eta < 0.800) return 0.0176;
  if (pt >= 20 && pt < 30 && eta >= 0.800 && eta < 1.442) return 0.0259;
  if (pt >= 20 && pt < 30 && eta >= 1.442 && eta < 1.566) return 0.1187;
  if (pt >= 20 && pt < 30 && eta >= 1.566 && eta < 2.000) return 0.0170;
  if (pt >= 20 && pt < 30 && eta >= 2.000 && eta < 2.500) return 0.0172;
  if (pt >= 30 && pt < 40 && eta >= -2.500 && eta < -2.000) return 0.0090;
  if (pt >= 30 && pt < 40 && eta >= -2.000 && eta < -1.566) return 0.0068;
  if (pt >= 30 && pt < 40 && eta >= -1.566 && eta < -1.442) return 0.0090;
  if (pt >= 30 && pt < 40 && eta >= -1.442 && eta < -0.800) return 0.0085;
  if (pt >= 30 && pt < 40 && eta >= -0.800 && eta < -0.000) return 0.0070;
  if (pt >= 30 && pt < 40 && eta >= -0.000 && eta < 0.800) return 0.0070;
  if (pt >= 30 && pt < 40 && eta >= 0.800 && eta < 1.442) return 0.0085;
  if (pt >= 30 && pt < 40 && eta >= 1.442 && eta < 1.566) return 0.0090;
  if (pt >= 30 && pt < 40 && eta >= 1.566 && eta < 2.000) return 0.0068;
  if (pt >= 30 && pt < 40 && eta >= 2.000 && eta < 2.500) return 0.0090;
  if (pt >= 40 && pt < 50 && eta >= -2.500 && eta < -2.000) return 0.0146;
  if (pt >= 40 && pt < 50 && eta >= -2.000 && eta < -1.566) return 0.0065;
  if (pt >= 40 && pt < 50 && eta >= -1.566 && eta < -1.442) return 0.0168;
  if (pt >= 40 && pt < 50 && eta >= -1.442 && eta < -0.800) return 0.0040;
  if (pt >= 40 && pt < 50 && eta >= -0.800 && eta < -0.000) return 0.0036;
  if (pt >= 40 && pt < 50 && eta >= -0.000 && eta < 0.800) return 0.0036;
  if (pt >= 40 && pt < 50 && eta >= 0.800 && eta < 1.442) return 0.0040;
  if (pt >= 40 && pt < 50 && eta >= 1.442 && eta < 1.566) return 0.0168;
  if (pt >= 40 && pt < 50 && eta >= 1.566 && eta < 2.000) return 0.0065;
  if (pt >= 40 && pt < 50 && eta >= 2.000 && eta < 2.500) return 0.0146;
  if (pt >= 50 && pt < 100 && eta >= -2.500 && eta < -2.000) return 0.0195;
  if (pt >= 50 && pt < 100 && eta >= -2.000 && eta < -1.566) return 0.0111;
  if (pt >= 50 && pt < 100 && eta >= -1.566 && eta < -1.442) return 0.0171;
  if (pt >= 50 && pt < 100 && eta >= -1.442 && eta < -0.800) return 0.0057;
  if (pt >= 50 && pt < 100 && eta >= -0.800 && eta < -0.000) return 0.0085;
  if (pt >= 50 && pt < 100 && eta >= -0.000 && eta < 0.800) return 0.0085;
  if (pt >= 50 && pt < 100 && eta >= 0.800 && eta < 1.442) return 0.0057;
  if (pt >= 50 && pt < 100 && eta >= 1.442 && eta < 1.566) return 0.0171;
  if (pt >= 50 && pt < 100 && eta >= 1.566 && eta < 2.000) return 0.0111;
  if (pt >= 50 && pt < 100 && eta >= 2.000 && eta < 2.500) return 0.0195;
  if (pt >= 100  && eta >= -2.500 && eta < -2.000) return 0.0417;
  if (pt >= 100  && eta >= -2.000 && eta < -1.566) return 0.0146;
  if (pt >= 100  && eta >= -1.566 && eta < -1.442) return 0.0636;
  if (pt >= 100  && eta >= -1.442 && eta < -0.800) return 0.0084;
  if (pt >= 100  && eta >= -0.800 && eta < -0.000) return 0.0062;
  if (pt >= 100  && eta >= -0.000 && eta < 0.800) return 0.0062;
  if (pt >= 100  && eta >= 0.800 && eta < 1.442) return 0.0084;
  if (pt >= 100  && eta >= 1.442 && eta < 1.566) return 0.0636;
  if (pt >= 100  && eta >= 1.566 && eta < 2.000) return 0.0146;
  if (pt >= 100  && eta >= 2.000 && eta < 2.500) return 0.0417;
  return 1.;
}


float electronScaleFactor_legacy(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 1.0634;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.9313;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.444) return 1.0123;
  if (pt >= 10 && pt < 20 && eta >= -1.444 && eta < -0.800) return 0.9516;
  if (pt >= 10 && pt < 20 && eta >= -0.800 && eta < 0.000) return 0.9588;
  if (pt >= 10 && pt < 20 && eta >= 0.000 && eta < 0.800) return 0.9735;
  if (pt >= 10 && pt < 20 && eta >= 0.800 && eta < 1.444) return 0.9210;
  if (pt >= 10 && pt < 20 && eta >= 1.444 && eta < 1.566) return 1.0384;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.9239;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 0.9779;
  if (pt >= 20 && pt < 35 && eta >= -2.500 && eta < -2.000) return 0.9979;
  if (pt >= 20 && pt < 35 && eta >= -2.000 && eta < -1.566) return 0.9300;
  if (pt >= 20 && pt < 35 && eta >= -1.566 && eta < -1.444) return 0.9336;
  if (pt >= 20 && pt < 35 && eta >= -1.444 && eta < -0.800) return 0.9176;
  if (pt >= 20 && pt < 35 && eta >= -0.800 && eta < 0.000) return 0.9277;
  if (pt >= 20 && pt < 35 && eta >= 0.000 && eta < 0.800) return 0.9426;
  if (pt >= 20 && pt < 35 && eta >= 0.800 && eta < 1.444) return 0.8963;
  if (pt >= 20 && pt < 35 && eta >= 1.444 && eta < 1.566) return 0.8672;
  if (pt >= 20 && pt < 35 && eta >= 1.566 && eta < 2.000) return 0.9107;
  if (pt >= 20 && pt < 35 && eta >= 2.000 && eta < 2.500) return 0.9431;
  if (pt >= 35 && pt < 50 && eta >= -2.500 && eta < -2.000) return 0.9861;
  if (pt >= 35 && pt < 50 && eta >= -2.000 && eta < -1.566) return 0.9431;
  if (pt >= 35 && pt < 50 && eta >= -1.566 && eta < -1.444) return 0.9259;
  if (pt >= 35 && pt < 50 && eta >= -1.444 && eta < -0.800) return 0.9339;
  if (pt >= 35 && pt < 50 && eta >= -0.800 && eta < 0.000) return 0.9336;
  if (pt >= 35 && pt < 50 && eta >= 0.000 && eta < 0.800) return 0.9471;
  if (pt >= 35 && pt < 50 && eta >= 0.800 && eta < 1.444) return 0.9155;
  if (pt >= 35 && pt < 50 && eta >= 1.444 && eta < 1.566) return 0.9029;
  if (pt >= 35 && pt < 50 && eta >= 1.566 && eta < 2.000) return 0.9310;
  if (pt >= 35 && pt < 50 && eta >= 2.000 && eta < 2.500) return 0.9554;
  if (pt >= 50 && pt < 100 && eta >= -2.500 && eta < -2.000) return 0.9949;
  if (pt >= 50 && pt < 100 && eta >= -2.000 && eta < -1.566) return 0.9619;
  if (pt >= 50 && pt < 100 && eta >= -1.566 && eta < -1.444) return 0.9366;
  if (pt >= 50 && pt < 100 && eta >= -1.444 && eta < -0.800) return 0.9334;
  if (pt >= 50 && pt < 100 && eta >= -0.800 && eta < 0.000) return 0.9351;
  if (pt >= 50 && pt < 100 && eta >= 0.000 && eta < 0.800) return 0.9480;
  if (pt >= 50 && pt < 100 && eta >= 0.800 && eta < 1.444) return 0.9185;
  if (pt >= 50 && pt < 100 && eta >= 1.444 && eta < 1.566) return 0.9167;
  if (pt >= 50 && pt < 100 && eta >= 1.566 && eta < 2.000) return 0.9569;
  if (pt >= 50 && pt < 100 && eta >= 2.000 && eta < 2.500) return 0.9662;
  if (pt >= 100 && pt < 200 && eta >= -2.500 && eta < -2.000) return 1.0214;
  if (pt >= 100 && pt < 200 && eta >= -2.000 && eta < -1.566) return 0.9882;
  if (pt >= 100 && pt < 200 && eta >= -1.566 && eta < -1.444) return 0.9216;
  if (pt >= 100 && pt < 200 && eta >= -1.444 && eta < -0.800) return 0.9215;
  if (pt >= 100 && pt < 200 && eta >= -0.800 && eta < 0.000) return 0.9382;
  if (pt >= 100 && pt < 200 && eta >= 0.000 && eta < 0.800) return 0.9537;
  if (pt >= 100 && pt < 200 && eta >= 0.800 && eta < 1.444) return 0.9074;
  if (pt >= 100 && pt < 200 && eta >= 1.444 && eta < 1.566) return 0.8905;
  if (pt >= 100 && pt < 200 && eta >= 1.566 && eta < 2.000) return 0.9314;
  if (pt >= 100 && pt < 200 && eta >= 2.000 && eta < 2.500) return 1.0301;
  if (pt >= 200  && eta >= -2.500 && eta < -2.000) return 0.7836;
  if (pt >= 200  && eta >= -2.000 && eta < -1.566) return 0.9927;
  if (pt >= 200  && eta >= -1.566 && eta < -1.444) return 1.0170;
  if (pt >= 200  && eta >= -1.444 && eta < -0.800) return 0.9272;
  if (pt >= 200  && eta >= -0.800 && eta < 0.000) return 0.9366;
  if (pt >= 200  && eta >= 0.000 && eta < 0.800) return 0.9238;
  if (pt >= 200  && eta >= 0.800 && eta < 1.444) return 0.9028;
  if (pt >= 200  && eta >= 1.444 && eta < 1.566) return 0.8326;
  if (pt >= 200  && eta >= 1.566 && eta < 2.000) return 0.9973;
  if (pt >= 200  && eta >= 2.000 && eta < 2.500) return 1.0125;
  return 0.0;
}

float electronScaleFactorError_legacy(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 0.0366;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.0226;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.444) return 0.0666;
  if (pt >= 10 && pt < 20 && eta >= -1.444 && eta < -0.800) return 0.0215;
  if (pt >= 10 && pt < 20 && eta >= -0.800 && eta < 0.000) return 0.0233;
  if (pt >= 10 && pt < 20 && eta >= 0.000 && eta < 0.800) return 0.0236;
  if (pt >= 10 && pt < 20 && eta >= 0.800 && eta < 1.444) return 0.0214;
  if (pt >= 10 && pt < 20 && eta >= 1.444 && eta < 1.566) return 0.0663;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.0226;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 0.0350;
  if (pt >= 20 && pt < 35 && eta >= -2.500 && eta < -2.000) return 0.0201;
  if (pt >= 20 && pt < 35 && eta >= -2.000 && eta < -1.566) return 0.0192;
  if (pt >= 20 && pt < 35 && eta >= -1.566 && eta < -1.444) return 0.0681;
  if (pt >= 20 && pt < 35 && eta >= -1.444 && eta < -0.800) return 0.0219;
  if (pt >= 20 && pt < 35 && eta >= -0.800 && eta < 0.000) return 0.0155;
  if (pt >= 20 && pt < 35 && eta >= 0.000 && eta < 0.800) return 0.0156;
  if (pt >= 20 && pt < 35 && eta >= 0.800 && eta < 1.444) return 0.0216;
  if (pt >= 20 && pt < 35 && eta >= 1.444 && eta < 1.566) return 0.0661;
  if (pt >= 20 && pt < 35 && eta >= 1.566 && eta < 2.000) return 0.0190;
  if (pt >= 20 && pt < 35 && eta >= 2.000 && eta < 2.500) return 0.0194;
  if (pt >= 35 && pt < 50 && eta >= -2.500 && eta < -2.000) return 0.0085;
  if (pt >= 35 && pt < 50 && eta >= -2.000 && eta < -1.566) return 0.0069;
  if (pt >= 35 && pt < 50 && eta >= -1.566 && eta < -1.444) return 0.0115;
  if (pt >= 35 && pt < 50 && eta >= -1.444 && eta < -0.800) return 0.0062;
  if (pt >= 35 && pt < 50 && eta >= -0.800 && eta < 0.000) return 0.0047;
  if (pt >= 35 && pt < 50 && eta >= 0.000 && eta < 0.800) return 0.0048;
  if (pt >= 35 && pt < 50 && eta >= 0.800 && eta < 1.444) return 0.0062;
  if (pt >= 35 && pt < 50 && eta >= 1.444 && eta < 1.566) return 0.0114;
  if (pt >= 35 && pt < 50 && eta >= 1.566 && eta < 2.000) return 0.0069;
  if (pt >= 35 && pt < 50 && eta >= 2.000 && eta < 2.500) return 0.0084;
  if (pt >= 50 && pt < 100 && eta >= -2.500 && eta < -2.000) return 0.0282;
  if (pt >= 50 && pt < 100 && eta >= -2.000 && eta < -1.566) return 0.0135;
  if (pt >= 50 && pt < 100 && eta >= -1.566 && eta < -1.444) return 0.0224;
  if (pt >= 50 && pt < 100 && eta >= -1.444 && eta < -0.800) return 0.0132;
  if (pt >= 50 && pt < 100 && eta >= -0.800 && eta < 0.000) return 0.0140;
  if (pt >= 50 && pt < 100 && eta >= 0.000 && eta < 0.800) return 0.0140;
  if (pt >= 50 && pt < 100 && eta >= 0.800 && eta < 1.444) return 0.0131;
  if (pt >= 50 && pt < 100 && eta >= 1.444 && eta < 1.566) return 0.0222;
  if (pt >= 50 && pt < 100 && eta >= 1.566 && eta < 2.000) return 0.0135;
  if (pt >= 50 && pt < 100 && eta >= 2.000 && eta < 2.500) return 0.0275;
  if (pt >= 100 && pt < 200 && eta >= -2.500 && eta < -2.000) return 0.0383;
  if (pt >= 100 && pt < 200 && eta >= -2.000 && eta < -1.566) return 0.0237;
  if (pt >= 100 && pt < 200 && eta >= -1.566 && eta < -1.444) return 0.0431;
  if (pt >= 100 && pt < 200 && eta >= -1.444 && eta < -0.800) return 0.0143;
  if (pt >= 100 && pt < 200 && eta >= -0.800 && eta < 0.000) return 0.0073;
  if (pt >= 100 && pt < 200 && eta >= 0.000 && eta < 0.800) return 0.0075;
  if (pt >= 100 && pt < 200 && eta >= 0.800 && eta < 1.444) return 0.0141;
  if (pt >= 100 && pt < 200 && eta >= 1.444 && eta < 1.566) return 0.0420;
  if (pt >= 100 && pt < 200 && eta >= 1.566 && eta < 2.000) return 0.0230;
  if (pt >= 100 && pt < 200 && eta >= 2.000 && eta < 2.500) return 0.0386;
  if (pt >= 200  && eta >= -2.500 && eta < -2.000) return 0.1416;
  if (pt >= 200  && eta >= -2.000 && eta < -1.566) return 0.0688;
  if (pt >= 200  && eta >= -1.566 && eta < -1.444) return 0.1483;
  if (pt >= 200  && eta >= -1.444 && eta < -0.800) return 0.0275;
  if (pt >= 200  && eta >= -0.800 && eta < 0.000) return 0.0352;
  if (pt >= 200  && eta >= 0.000 && eta < 0.800) return 0.0347;
  if (pt >= 200  && eta >= 0.800 && eta < 1.444) return 0.0269;
  if (pt >= 200  && eta >= 1.444 && eta < 1.566) return 0.1274;
  if (pt >= 200  && eta >= 1.566 && eta < 2.000) return 0.0703;
  if (pt >= 200  && eta >= 2.000 && eta < 2.500) return 0.1623;
  return 0.0;
}

float electronScaleFactorReco_legacy(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 1.0423;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.9744;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.444) return 1.4277;
  if (pt >= 10 && pt < 20 && eta >= -1.444 && eta < -1.000) return 0.9894;
  if (pt >= 10 && pt < 20 && eta >= -1.000 && eta < 0.000) return 0.9820;
  if (pt >= 10 && pt < 20 && eta >= 0.000 && eta < 1.000) return 0.9820;
  if (pt >= 10 && pt < 20 && eta >= 1.000 && eta < 1.444) return 0.9894;
  if (pt >= 10 && pt < 20 && eta >= 1.444 && eta < 1.566) return 1.4277;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.9744;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 1.0423;
  if (pt >= 20 && pt < 45 && eta >= -2.500 && eta < -2.000) return 1.0164;
  if (pt >= 20 && pt < 45 && eta >= -2.000 && eta < -1.566) return 0.9979;
  if (pt >= 20 && pt < 45 && eta >= -1.566 && eta < -1.444) return 0.9908;
  if (pt >= 20 && pt < 45 && eta >= -1.444 && eta < -1.000) return 0.9918;
  if (pt >= 20 && pt < 45 && eta >= -1.000 && eta < -0.500) return 0.9867;
  if (pt >= 20 && pt < 45 && eta >= -0.500 && eta < 0.000) return 0.9836;
  if (pt >= 20 && pt < 45 && eta >= 0.000 && eta < 0.500) return 0.9836;
  if (pt >= 20 && pt < 45 && eta >= 0.500 && eta < 1.000) return 0.9867;
  if (pt >= 20 && pt < 45 && eta >= 1.000 && eta < 1.444) return 0.9918;
  if (pt >= 20 && pt < 45 && eta >= 1.444 && eta < 1.566) return 0.9908;
  if (pt >= 20 && pt < 45 && eta >= 1.566 && eta < 2.000) return 0.9979;
  if (pt >= 20 && pt < 45 && eta >= 2.000 && eta < 2.500) return 1.0164;
  if (pt >= 45 && pt < 75 && eta >= -2.500 && eta < -2.000) return 1.0021;
  if (pt >= 45 && pt < 75 && eta >= -2.000 && eta < -1.566) return 0.9969;
  if (pt >= 45 && pt < 75 && eta >= -1.566 && eta < -1.444) return 0.9622;
  if (pt >= 45 && pt < 75 && eta >= -1.444 && eta < -1.000) return 0.9918;
  if (pt >= 45 && pt < 75 && eta >= -1.000 && eta < -0.500) return 0.9878;
  if (pt >= 45 && pt < 75 && eta >= -0.500 && eta < 0.000) return 0.9868;
  if (pt >= 45 && pt < 75 && eta >= 0.000 && eta < 0.500) return 0.9868;
  if (pt >= 45 && pt < 75 && eta >= 0.500 && eta < 1.000) return 0.9878;
  if (pt >= 45 && pt < 75 && eta >= 1.000 && eta < 1.444) return 0.9918;
  if (pt >= 45 && pt < 75 && eta >= 1.444 && eta < 1.566) return 0.9622;
  if (pt >= 45 && pt < 75 && eta >= 1.566 && eta < 2.000) return 0.9969;
  if (pt >= 45 && pt < 75 && eta >= 2.000 && eta < 2.500) return 1.0021;
  if (pt >= 75 && pt < 100 && eta >= -2.500 && eta < -2.000) return 1.0180;
  if (pt >= 75 && pt < 100 && eta >= -2.000 && eta < -1.566) return 1.0154;
  if (pt >= 75 && pt < 100 && eta >= -1.566 && eta < -1.444) return 1.0329;
  if (pt >= 75 && pt < 100 && eta >= -1.444 && eta < -1.000) return 1.0081;
  if (pt >= 75 && pt < 100 && eta >= -1.000 && eta < -0.500) return 1.0051;
  if (pt >= 75 && pt < 100 && eta >= -0.500 && eta < 0.000) return 0.9969;
  if (pt >= 75 && pt < 100 && eta >= 0.000 && eta < 0.500) return 0.9969;
  if (pt >= 75 && pt < 100 && eta >= 0.500 && eta < 1.000) return 1.0051;
  if (pt >= 75 && pt < 100 && eta >= 1.000 && eta < 1.444) return 1.0081;
  if (pt >= 75 && pt < 100 && eta >= 1.444 && eta < 1.566) return 1.0329;
  if (pt >= 75 && pt < 100 && eta >= 1.566 && eta < 2.000) return 1.0154;
  if (pt >= 75 && pt < 100 && eta >= 2.000 && eta < 2.500) return 1.0180;
  if (pt >= 100  && eta >= -2.500 && eta < -2.000) return 0.9843;
  if (pt >= 100  && eta >= -2.000 && eta < -1.566) return 1.0000;
  if (pt >= 100  && eta >= -1.566 && eta < -1.444) return 1.0022;
  if (pt >= 100  && eta >= -1.444 && eta < -1.000) return 0.9869;
  if (pt >= 100  && eta >= -1.000 && eta < -0.500) return 0.9939;
  if (pt >= 100  && eta >= -0.500 && eta < 0.000) return 0.9858;
  if (pt >= 100  && eta >= 0.000 && eta < 0.500) return 0.9858;
  if (pt >= 100  && eta >= 0.500 && eta < 1.000) return 0.9939;
  if (pt >= 100  && eta >= 1.000 && eta < 1.444) return 0.9869;
  if (pt >= 100  && eta >= 1.444 && eta < 1.566) return 1.0022;
  if (pt >= 100  && eta >= 1.566 && eta < 2.000) return 1.0000;
  if (pt >= 100  && eta >= 2.000 && eta < 2.500) return 0.9843;
  return 0.0;
}

float electronScaleFactorRecoError_legacy(float pt, float eta) {
  if (pt >= 10 && pt < 20 && eta >= -2.500 && eta < -2.000) return 0.0142;
  if (pt >= 10 && pt < 20 && eta >= -2.000 && eta < -1.566) return 0.0103;
  if (pt >= 10 && pt < 20 && eta >= -1.566 && eta < -1.444) return 0.1164;
  if (pt >= 10 && pt < 20 && eta >= -1.444 && eta < -1.000) return 0.0190;
  if (pt >= 10 && pt < 20 && eta >= -1.000 && eta < 0.000) return 0.0226;
  if (pt >= 10 && pt < 20 && eta >= 0.000 && eta < 1.000) return 0.0226;
  if (pt >= 10 && pt < 20 && eta >= 1.000 && eta < 1.444) return 0.0190;
  if (pt >= 10 && pt < 20 && eta >= 1.444 && eta < 1.566) return 0.1164;
  if (pt >= 10 && pt < 20 && eta >= 1.566 && eta < 2.000) return 0.0103;
  if (pt >= 10 && pt < 20 && eta >= 2.000 && eta < 2.500) return 0.0142;
  if (pt >= 20 && pt < 45 && eta >= -2.500 && eta < -2.000) return 0.0039;
  if (pt >= 20 && pt < 45 && eta >= -2.000 && eta < -1.566) return 0.0027;
  if (pt >= 20 && pt < 45 && eta >= -1.566 && eta < -1.444) return 0.0799;
  if (pt >= 20 && pt < 45 && eta >= -1.444 && eta < -1.000) return 0.0027;
  if (pt >= 20 && pt < 45 && eta >= -1.000 && eta < -0.500) return 0.0025;
  if (pt >= 20 && pt < 45 && eta >= -0.500 && eta < 0.000) return 0.0027;
  if (pt >= 20 && pt < 45 && eta >= 0.000 && eta < 0.500) return 0.0027;
  if (pt >= 20 && pt < 45 && eta >= 0.500 && eta < 1.000) return 0.0025;
  if (pt >= 20 && pt < 45 && eta >= 1.000 && eta < 1.444) return 0.0027;
  if (pt >= 20 && pt < 45 && eta >= 1.444 && eta < 1.566) return 0.0799;
  if (pt >= 20 && pt < 45 && eta >= 1.566 && eta < 2.000) return 0.0027;
  if (pt >= 20 && pt < 45 && eta >= 2.000 && eta < 2.500) return 0.0039;
  if (pt >= 45 && pt < 75 && eta >= -2.500 && eta < -2.000) return 0.0064;
  if (pt >= 45 && pt < 75 && eta >= -2.000 && eta < -1.566) return 0.0045;
  if (pt >= 45 && pt < 75 && eta >= -1.566 && eta < -1.444) return 0.0230;
  if (pt >= 45 && pt < 75 && eta >= -1.444 && eta < -1.000) return 0.0053;
  if (pt >= 45 && pt < 75 && eta >= -1.000 && eta < -0.500) return 0.0035;
  if (pt >= 45 && pt < 75 && eta >= -0.500 && eta < 0.000) return 0.0035;
  if (pt >= 45 && pt < 75 && eta >= 0.000 && eta < 0.500) return 0.0035;
  if (pt >= 45 && pt < 75 && eta >= 0.500 && eta < 1.000) return 0.0035;
  if (pt >= 45 && pt < 75 && eta >= 1.000 && eta < 1.444) return 0.0053;
  if (pt >= 45 && pt < 75 && eta >= 1.444 && eta < 1.566) return 0.0230;
  if (pt >= 45 && pt < 75 && eta >= 1.566 && eta < 2.000) return 0.0045;
  if (pt >= 45 && pt < 75 && eta >= 2.000 && eta < 2.500) return 0.0064;
  if (pt >= 75 && pt < 100 && eta >= -2.500 && eta < -2.000) return 0.0122;
  if (pt >= 75 && pt < 100 && eta >= -2.000 && eta < -1.566) return 0.0107;
  if (pt >= 75 && pt < 100 && eta >= -1.566 && eta < -1.444) return 0.0441;
  if (pt >= 75 && pt < 100 && eta >= -1.444 && eta < -1.000) return 0.0056;
  if (pt >= 75 && pt < 100 && eta >= -1.000 && eta < -0.500) return 0.0047;
  if (pt >= 75 && pt < 100 && eta >= -0.500 && eta < 0.000) return 0.0064;
  if (pt >= 75 && pt < 100 && eta >= 0.000 && eta < 0.500) return 0.0064;
  if (pt >= 75 && pt < 100 && eta >= 0.500 && eta < 1.000) return 0.0047;
  if (pt >= 75 && pt < 100 && eta >= 1.000 && eta < 1.444) return 0.0056;
  if (pt >= 75 && pt < 100 && eta >= 1.444 && eta < 1.566) return 0.0441;
  if (pt >= 75 && pt < 100 && eta >= 1.566 && eta < 2.000) return 0.0107;
  if (pt >= 75 && pt < 100 && eta >= 2.000 && eta < 2.500) return 0.0122;
  if (pt >= 100  && eta >= -2.500 && eta < -2.000) return 0.0174;
  if (pt >= 100  && eta >= -2.000 && eta < -1.566) return 0.0102;
  if (pt >= 100  && eta >= -1.566 && eta < -1.444) return 0.0465;
  if (pt >= 100  && eta >= -1.444 && eta < -1.000) return 0.0068;
  if (pt >= 100  && eta >= -1.000 && eta < -0.500) return 0.0058;
  if (pt >= 100  && eta >= -0.500 && eta < 0.000) return 0.0069;
  if (pt >= 100  && eta >= 0.000 && eta < 0.500) return 0.0069;
  if (pt >= 100  && eta >= 0.500 && eta < 1.000) return 0.0058;
  if (pt >= 100  && eta >= 1.000 && eta < 1.444) return 0.0068;
  if (pt >= 100  && eta >= 1.444 && eta < 1.566) return 0.0465;
  if (pt >= 100  && eta >= 1.566 && eta < 2.000) return 0.0102;
  if (pt >= 100  && eta >= 2.000 && eta < 2.500) return 0.0174;
  return 0.0;
}

