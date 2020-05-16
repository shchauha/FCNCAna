//sample: 1-ttH
//sample: 2-tHW
//sample: 3-tHq
//sample += 3 for corresponding pseudoscalar
int nPoints_higgs(int sample, int mH){
    if (sample == 1){
        if (mH == 350) return 103778;
        if (mH == 370) return 91175;
        if (mH == 390) return 93998;
        if (mH == 410) return 98687;
        if (mH == 430) return 102421;
        if (mH == 450) return 87667;
        if (mH == 470) return 96639;
        if (mH == 490) return 93382;
        if (mH == 510) return 101248;
        if (mH == 530) return 86954;
        if (mH == 550) return 91983;
        if (mH == 600) return 97030;
        if (mH == 650) return 88140;
        if (mH == 700) return 92904;
        if (mH == 750) return 95115;
        if (mH == 800) return 86535;
        if (mH == 850) return 96488;
        if (mH == 900) return 93990;
    }
    else if (sample == 2){
        if (mH == 350) return 99550;
        if (mH == 370) return 100616;
        if (mH == 390) return 99218;
        if (mH == 410) return 104747;
        if (mH == 430) return 97215;
        if (mH == 450) return 100081;
        if (mH == 470) return 104684;
        if (mH == 490) return 90082;
        if (mH == 510) return 99372;
        if (mH == 530) return 103017;
        if (mH == 550) return 95217;
    }
    else if (sample == 3){
        if (mH == 350) return 95842;
        if (mH == 370) return 104916;
        if (mH == 390) return 97056;
        if (mH == 410) return 105446;
        if (mH == 430) return 100314;
        if (mH == 450) return 98848;
        if (mH == 470) return 95518;
        if (mH == 490) return 107314;
        if (mH == 510) return 98318;
        if (mH == 530) return 94914;
        if (mH == 550) return 93782;
    }
    return -99999999; 
}

////sample += 3 for corresponding pseudoscalar
//float xsec_higgs(int sample, int mH){
//    if (sample == 1){
//        if (mH == 350) return 0.001*24.2;
//        if (mH == 370) return 0.001*20.1;
//        if (mH == 390) return 0.001*17.8;
//        if (mH == 410) return 0.001*15.8;
//        if (mH == 430) return 0.001*14.9;
//        if (mH == 450) return 0.001*11.9;
//        if (mH == 470) return 0.001*11.2;
//        if (mH == 490) return 0.001*9.87;
//        if (mH == 510) return 0.001*9.75;
//        if (mH == 530) return 0.001*7.78;
//        if (mH == 550) return 0.001*7.49;
//        if (mH == 570) return 0.007207;
//        if (mH == 590) return 0.006467;
//        if (mH == 610) return 0.005803;
//        if (mH == 630) return 0.005222;
//        if (mH == 650) return 0.004728;
//        if (mH == 670) return 0.004325;
//        if (mH == 690) return 0.004020;
//        // if (mH == 600) return 0.00578;
//        // if (mH == 650) return 0.00465;
//        // if (mH == 700) return 0.00364;
//        // if (mH == 750) return 0.00299;
//        // if (mH == 800) return 0.00233;
//        // if (mH == 850) return 0.00195;
//        // if (mH == 900) return 0.00160;
//    }
//    else if (sample == 2){
//        if (mH == 350) return 0.001*8.37;
//        if (mH == 370) return 0.001*7.61;
//        if (mH == 390) return 0.001*7.07;
//        if (mH == 410) return 0.001*6.58;
//        if (mH == 430) return 0.001*5.72;
//        if (mH == 450) return 0.001*5.16;
//        if (mH == 470) return 0.001*4.81;
//        if (mH == 490) return 0.001*4.34;
//        if (mH == 510) return 0.001*4.19;
//        if (mH == 530) return 0.001*3.86;
//        if (mH == 550) return 0.001*3.41;
//        if (mH == 570) return 0.003378;
//        if (mH == 590) return 0.003095;
//        if (mH == 610) return 0.002830;
//        if (mH == 630) return 0.002585;
//        if (mH == 650) return 0.002361;
//        if (mH == 670) return 0.002157;
//        if (mH == 690) return 0.001976;
//    }
//    else if (sample == 3){
//        if (mH == 350) return 0.001*15.8;
//        if (mH == 370) return 0.001*13.5;
//        if (mH == 390) return 0.001*11.7;
//        if (mH == 410) return 0.001*10.2;
//        if (mH == 430) return 0.001*8.82;
//        if (mH == 450) return 0.001*7.73;
//        if (mH == 470) return 0.001*6.78;
//        if (mH == 490) return 0.001*5.93;
//        if (mH == 510) return 0.001*5.23;
//        if (mH == 530) return 0.001*4.59;
//        if (mH == 550) return 0.001*4.11;
//        if (mH == 570) return 0.003786;
//        if (mH == 590) return 0.003403;
//        if (mH == 610) return 0.003071;
//        if (mH == 630) return 0.002778;
//        if (mH == 650) return 0.002514;
//        if (mH == 670) return 0.002267;
//        if (mH == 690) return 0.002024;
//    }
//    if (sample == 4){
//        if (mH == 350) return 0.001*33.7;
//        if (mH == 370) return 0.001*28.6;
//        if (mH == 390) return 0.001*25.3;
//        if (mH == 410) return 0.001*22.0;
//        if (mH == 430) return 0.001*19.9;
//        if (mH == 450) return 0.001*17.7;
//        if (mH == 470) return 0.001*15.1;
//        if (mH == 490) return 0.001*14.9;
//        if (mH == 510) return 0.001*12.6;
//        if (mH == 530) return 0.001*11.5;
//        if (mH == 550) return 0.001*10.0;
//        if (mH == 570) return 0.009233;
//        if (mH == 590) return 0.008270;
//        if (mH == 610) return 0.007430;
//        if (mH == 630) return 0.006705;
//        if (mH == 650) return 0.006086;
//        if (mH == 670) return 0.005564;
//        if (mH == 690) return 0.005131;
//    }
//    else if (sample == 5){
//        if (mH == 350) return 0.001*7.97;
//        if (mH == 370) return 0.001*7.09;
//        if (mH == 390) return 0.001*6.88;
//        if (mH == 410) return 0.001*5.92;
//        if (mH == 430) return 0.001*5.39;
//        if (mH == 450) return 0.001*4.97;
//        if (mH == 470) return 0.001*4.59;
//        if (mH == 490) return 0.001*4.27;
//        if (mH == 510) return 0.001*3.87;
//        if (mH == 530) return 0.001*3.82;
//        if (mH == 550) return 0.001*3.26;
//        if (mH == 570) return 0.003167;
//        if (mH == 590) return 0.002946;
//        if (mH == 610) return 0.002739;
//        if (mH == 630) return 0.002542;
//        if (mH == 650) return 0.002348;
//        if (mH == 670) return 0.002155;
//        if (mH == 690) return 0.001958;
//    }
//    else if (sample == 6){
//        if (mH == 350) return 0.001*9.28;
//        if (mH == 370) return 0.001*8.14;
//        if (mH == 390) return 0.001*7.19;
//        if (mH == 410) return 0.001*6.43;
//        if (mH == 430) return 0.001*5.72;
//        if (mH == 450) return 0.001*5.07;
//        if (mH == 470) return 0.001*4.55;
//        if (mH == 490) return 0.001*4.09;
//        if (mH == 510) return 0.001*3.67;
//        if (mH == 530) return 0.001*3.27;
//        if (mH == 550) return 0.001*2.94;
//        if (mH == 570) return 0.002753;
//        if (mH == 590) return 0.002508;
//        if (mH == 610) return 0.002289;
//        if (mH == 630) return 0.002090;
//        if (mH == 650) return 0.001903;
//        if (mH == 670) return 0.001723;
//        if (mH == 690) return 0.001543;
//    }
//    return -1.0;
//}


// sample: 1-ttH,2-tHW,3-tHq; +3 for corresponding pseudoscalar
float xsec_higgs(int sample, int mH) {
    if (sample == 1) {
        if (mH == 350) return 0.01851;
        if (mH == 370) return 0.01698;
        if (mH == 390) return 0.01510;
        if (mH == 410) return 0.01343;
        if (mH == 430) return 0.01194;
        if (mH == 450) return 0.01068;
        if (mH == 470) return 0.00961;
        if (mH == 490) return 0.00863;
        if (mH == 510) return 0.00777;
        if (mH == 530) return 0.00700;
        if (mH == 550) return 0.00634;
        if (mH == 570) return 0.00575;
        if (mH == 590) return 0.00520;
        if (mH == 610) return 0.00473;
        if (mH == 630) return 0.00429;
        if (mH == 650) return 0.00391;
    }
    else if (sample == 2) {
        if (mH == 350) return 0.00866;
        if (mH == 370) return 0.00819;
        if (mH == 390) return 0.00749;
        if (mH == 410) return 0.00685;
        if (mH == 430) return 0.00624;
        if (mH == 450) return 0.00568;
        if (mH == 470) return 0.00518;
        if (mH == 490) return 0.00470;
        if (mH == 510) return 0.00428;
        if (mH == 530) return 0.00391;
        if (mH == 550) return 0.00356;
        if (mH == 570) return 0.00325;
        if (mH == 590) return 0.00298;
        if (mH == 610) return 0.00272;
        if (mH == 630) return 0.00250;
        if (mH == 650) return 0.00228;
    }
    else if (sample == 3) {
        if (mH == 350) return 0.01774;
        if (mH == 370) return 0.01588;
        if (mH == 390) return 0.01389;
        if (mH == 410) return 0.01217;
        if (mH == 430) return 0.01063;
        if (mH == 450) return 0.00932;
        if (mH == 470) return 0.00822;
        if (mH == 490) return 0.00719;
        if (mH == 510) return 0.00635;
        if (mH == 530) return 0.00562;
        if (mH == 550) return 0.00500;
        if (mH == 570) return 0.00444;
        if (mH == 590) return 0.00396;
        if (mH == 610) return 0.00354;
        if (mH == 630) return 0.00316;
        if (mH == 650) return 0.00283;
    }
    else if (sample == 4) {
        if (mH == 350) return 0.02899;
        if (mH == 370) return 0.02564;
        if (mH == 390) return 0.02262;
        if (mH == 410) return 0.01990;
        if (mH == 430) return 0.01759;
        if (mH == 450) return 0.01552;
        if (mH == 470) return 0.01378;
        if (mH == 490) return 0.01225;
        if (mH == 510) return 0.01093;
        if (mH == 530) return 0.00974;
        if (mH == 550) return 0.00871;
        if (mH == 570) return 0.00780;
        if (mH == 590) return 0.00699;
        if (mH == 610) return 0.00627;
        if (mH == 630) return 0.00563;
        if (mH == 650) return 0.00508;
    }
    else if (sample == 5) {
        if (mH == 350) return 0.00898;
        if (mH == 370) return 0.00818;
        if (mH == 390) return 0.00743;
        if (mH == 410) return 0.00676;
        if (mH == 430) return 0.00613;
        if (mH == 450) return 0.00556;
        if (mH == 470) return 0.00506;
        if (mH == 490) return 0.00461;
        if (mH == 510) return 0.00420;
        if (mH == 530) return 0.00385;
        if (mH == 550) return 0.00351;
        if (mH == 570) return 0.00321;
        if (mH == 590) return 0.00294;
        if (mH == 610) return 0.00269;
        if (mH == 630) return 0.00247;
        if (mH == 650) return 0.00225;
    }
    else if (sample == 6) {
        if (mH == 350) return 0.01130;
        if (mH == 370) return 0.01002;
        if (mH == 390) return 0.00885;
        if (mH == 410) return 0.00787;
        if (mH == 430) return 0.00700;
        if (mH == 450) return 0.00626;
        if (mH == 470) return 0.00557;
        if (mH == 490) return 0.00498;
        if (mH == 510) return 0.00447;
        if (mH == 530) return 0.00401;
        if (mH == 550) return 0.00360;
        if (mH == 570) return 0.00323;
        if (mH == 590) return 0.00292;
        if (mH == 610) return 0.00263;
        if (mH == 630) return 0.00238;
        if (mH == 650) return 0.00215;
    }
    return -1;
}
