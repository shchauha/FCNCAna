double Ele23_Ele12_LeadingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >5&& pt<=6&& abseta>1.479&& abseta<=2) return 0;
  if(pt >5&& pt<=6&& abseta>2&& abseta<=2.5) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >6&& pt<=7&& abseta>1.479&& abseta<=2) return 0;
  if(pt >6&& pt<=7&& abseta>2&& abseta<=2.5) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.00028177;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.479) return 0.000220167;
  if(pt >7&& pt<=10&& abseta>1.479&& abseta<=2) return 0;
  if(pt >7&& pt<=10&& abseta>2&& abseta<=2.5) return 0;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.00431652;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.479) return 0.00291642;
  if(pt >10&& pt<=15&& abseta>1.479&& abseta<=2) return 0.00236546;
  if(pt >10&& pt<=15&& abseta>2&& abseta<=2.5) return 0.00170711;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.00542683;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.479) return 0.00659796;
  if(pt >15&& pt<=20&& abseta>1.479&& abseta<=2) return 0.0148951;
  if(pt >15&& pt<=20&& abseta>2&& abseta<=2.5) return 0.00375244;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.449631;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.479) return 0.391548;
  if(pt >20&& pt<=25&& abseta>1.479&& abseta<=2) return 0.463316;
  if(pt >20&& pt<=25&& abseta>2&& abseta<=2.5) return 0.319144;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.936224;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.479) return 0.945047;
  if(pt >25&& pt<=30&& abseta>1.479&& abseta<=2) return 0.971591;
  if(pt >25&& pt<=30&& abseta>2&& abseta<=2.5) return 0.931077;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.949822;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.479) return 0.963277;
  if(pt >30&& pt<=35&& abseta>1.479&& abseta<=2) return 0.981658;
  if(pt >30&& pt<=35&& abseta>2&& abseta<=2.5) return 0.967216;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.960488;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.479) return 0.969078;
  if(pt >35&& pt<=40&& abseta>1.479&& abseta<=2) return 0.985574;
  if(pt >35&& pt<=40&& abseta>2&& abseta<=2.5) return 0.976847;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.964928;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.479) return 0.972645;
  if(pt >40&& pt<=50&& abseta>1.479&& abseta<=2) return 0.988858;
  if(pt >40&& pt<=50&& abseta>2&& abseta<=2.5) return 0.981689;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.970105;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.479) return 0.974872;
  if(pt >50&& pt<=75&& abseta>1.479&& abseta<=2) return 0.988842;
  if(pt >50&& pt<=75&& abseta>2&& abseta<=2.5) return 0.986394;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.975481;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.479) return 0.980721;
  if(pt >75&& pt<=100&& abseta>1.479&& abseta<=2) return 0.993264;
  if(pt >75&& pt<=100&& abseta>2&& abseta<=2.5) return 0.992702;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.977292;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.479) return 0.981523;
  if(pt >100&& pt<=200&& abseta>1.479&& abseta<=2) return 0.99487;
  if(pt >100&& pt<=200&& abseta>2&& abseta<=2.5) return 0.99252;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.976461;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.479) return 0.984481;
  if(pt >200&& pt<=300&& abseta>1.479&& abseta<=2) return 0.996604;
  if(pt >200&& pt<=300&& abseta>2&& abseta<=2.5) return 0.997481;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.970077;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.479) return 0.983607;
  if(pt >300&& pt<=500&& abseta>1.479&& abseta<=2) return 0.99537;
  if(pt >300&& pt<=500&& abseta>2&& abseta<=2.5) return 1;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.970077;
  if(pt >500&& abseta>0.8&& abseta<=1.479) return 0.983607;
  if(pt >500&& abseta>1.479&& abseta<=2) return 0.99537;
  if(pt >500&& abseta>2&& abseta<=2.5) return 1;
  return 0;



}

double Ele23_Ele12_TrailingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >5&& pt<=6&& abseta>1.479&& abseta<=2) return 0;
  if(pt >5&& pt<=6&& abseta>2&& abseta<=2.5) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >6&& pt<=7&& abseta>1.479&& abseta<=2) return 0;
  if(pt >6&& pt<=7&& abseta>2&& abseta<=2.5) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.000422654;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.479) return 0.000440335;
  if(pt >7&& pt<=10&& abseta>1.479&& abseta<=2) return 0.000562272;
  if(pt >7&& pt<=10&& abseta>2&& abseta<=2.5) return 0;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.607299;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.479) return 0.568037;
  if(pt >10&& pt<=15&& abseta>1.479&& abseta<=2) return 0.500821;
  if(pt >10&& pt<=15&& abseta>2&& abseta<=2.5) return 0.431411;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.922124;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.479) return 0.937801;
  if(pt >15&& pt<=20&& abseta>1.479&& abseta<=2) return 0.925633;
  if(pt >15&& pt<=20&& abseta>2&& abseta<=2.5) return 0.843666;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.934702;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.479) return 0.947745;
  if(pt >20&& pt<=25&& abseta>1.479&& abseta<=2) return 0.97187;
  if(pt >20&& pt<=25&& abseta>2&& abseta<=2.5) return 0.908927;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.939428;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.479) return 0.955316;
  if(pt >25&& pt<=30&& abseta>1.479&& abseta<=2) return 0.980632;
  if(pt >25&& pt<=30&& abseta>2&& abseta<=2.5) return 0.946155;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.949931;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.479) return 0.963982;
  if(pt >30&& pt<=35&& abseta>1.479&& abseta<=2) return 0.982345;
  if(pt >30&& pt<=35&& abseta>2&& abseta<=2.5) return 0.967586;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.960491;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.479) return 0.969186;
  if(pt >35&& pt<=40&& abseta>1.479&& abseta<=2) return 0.985694;
  if(pt >35&& pt<=40&& abseta>2&& abseta<=2.5) return 0.977009;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.964929;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.479) return 0.972661;
  if(pt >40&& pt<=50&& abseta>1.479&& abseta<=2) return 0.988878;
  if(pt >40&& pt<=50&& abseta>2&& abseta<=2.5) return 0.981786;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.970105;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.479) return 0.974877;
  if(pt >50&& pt<=75&& abseta>1.479&& abseta<=2) return 0.988855;
  if(pt >50&& pt<=75&& abseta>2&& abseta<=2.5) return 0.986431;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.975481;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.479) return 0.980721;
  if(pt >75&& pt<=100&& abseta>1.479&& abseta<=2) return 0.993264;
  if(pt >75&& pt<=100&& abseta>2&& abseta<=2.5) return 0.992702;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.977292;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.479) return 0.981523;
  if(pt >100&& pt<=200&& abseta>1.479&& abseta<=2) return 0.99487;
  if(pt >100&& pt<=200&& abseta>2&& abseta<=2.5) return 0.99252;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.976461;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.479) return 0.984481;
  if(pt >200&& pt<=300&& abseta>1.479&& abseta<=2) return 0.996604;
  if(pt >200&& pt<=300&& abseta>2&& abseta<=2.5) return 0.997481;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.970077;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.479) return 0.983607;
  if(pt >300&& pt<=500&& abseta>1.479&& abseta<=2) return 0.99537;
  if(pt >300&& pt<=500&& abseta>2&& abseta<=2.5) return 1;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.970077;
  if(pt >500&& abseta>0.8&& abseta<=1.479) return 0.983607;
  if(pt >500&& abseta>1.479&& abseta<=2) return 0.99537;
  if(pt >500&& abseta>2&& abseta<=2.5) return 1;
  return 0;

  
}




double Ele23_Ele12_DZLeg_Eff(double eta1, double eta2){

  if(eta1 < eta2) {
    double buf =eta1;
    eta1=eta2;
    eta2=buf;
  }

  if(abs(eta1)>2.5) return 0.;
  if(abs(eta2)>2.5) return 0.;
  if(eta1 >-2.5&& eta1<=-2&& eta2>-2.5&& eta2<=-2) return 0.972408;
  if(eta1 >-2&& eta1<=-1.5&& eta2>-2.5&& eta2<=-2) return 0.976071;
  if(eta1 >-2&& eta1<=-1.5&& eta2>-2&& eta2<=-1.5) return 0.979875;
  if(eta1 >-1.5&& eta1<=-1&& eta2>-2.5&& eta2<=-2) return 0.982069;
  if(eta1 >-1.5&& eta1<=-1&& eta2>-2&& eta2<=-1.5) return 0.986621;
  if(eta1 >-1.5&& eta1<=-1&& eta2>-1.5&& eta2<=-1) return 0.99242;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-2.5&& eta2<=-2) return 0.983559;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-2&& eta2<=-1.5) return 0.987892;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-1.5&& eta2<=-1) return 0.994891;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-1&& eta2<=-0.5) return 0.996419;
  if(eta1 >-0.5&& eta1<=0&& eta2>-2.5&& eta2<=-2) return 0.981451;
  if(eta1 >-0.5&& eta1<=0&& eta2>-2&& eta2<=-1.5) return 0.987607;
  if(eta1 >-0.5&& eta1<=0&& eta2>-1.5&& eta2<=-1) return 0.994967;
  if(eta1 >-0.5&& eta1<=0&& eta2>-1&& eta2<=-0.5) return 0.996985;
  if(eta1 >-0.5&& eta1<=0&& eta2>-0.5&& eta2<=0) return 0.997191;
  if(eta1 >0&& eta1<=0.5&& eta2>-2.5&& eta2<=-2) return 0.984523;
  if(eta1 >0&& eta1<=0.5&& eta2>-2&& eta2<=-1.5) return 0.987602;
  if(eta1 >0&& eta1<=0.5&& eta2>-1.5&& eta2<=-1) return 0.994738;
  if(eta1 >0&& eta1<=0.5&& eta2>-1&& eta2<=-0.5) return 0.996844;
  if(eta1 >0&& eta1<=0.5&& eta2>-0.5&& eta2<=0) return 0.997149;
  if(eta1 >0&& eta1<=0.5&& eta2>0&& eta2<=0.5) return 0.997427;
  if(eta1 >0.5&& eta1<=1&& eta2>-2.5&& eta2<=-2) return 0.988543;
  if(eta1 >0.5&& eta1<=1&& eta2>-2&& eta2<=-1.5) return 0.988847;
  if(eta1 >0.5&& eta1<=1&& eta2>-1.5&& eta2<=-1) return 0.99475;
  if(eta1 >0.5&& eta1<=1&& eta2>-1&& eta2<=-0.5) return 0.997263;
  if(eta1 >0.5&& eta1<=1&& eta2>-0.5&& eta2<=0) return 0.997527;
  if(eta1 >0.5&& eta1<=1&& eta2>0&& eta2<=0.5) return 0.997728;
  if(eta1 >0.5&& eta1<=1&& eta2>0.5&& eta2<=1) return 0.998355;
  if(eta1 >1&& eta1<=1.5&& eta2>-2.5&& eta2<=-2) return 0.963636;
  if(eta1 >1&& eta1<=1.5&& eta2>-2&& eta2<=-1.5) return 0.987676;
  if(eta1 >1&& eta1<=1.5&& eta2>-1.5&& eta2<=-1) return 0.993714;
  if(eta1 >1&& eta1<=1.5&& eta2>-1&& eta2<=-0.5) return 0.995661;
  if(eta1 >1&& eta1<=1.5&& eta2>-0.5&& eta2<=0) return 0.995514;
  if(eta1 >1&& eta1<=1.5&& eta2>0&& eta2<=0.5) return 0.996184;
  if(eta1 >1&& eta1<=1.5&& eta2>0.5&& eta2<=1) return 0.996245;
  if(eta1 >1&& eta1<=1.5&& eta2>1&& eta2<=1.5) return 0.994036;
  if(eta1 >1.5&& eta1<=2&& eta2>-2.5&& eta2<=-2) return 0.972222;
  if(eta1 >1.5&& eta1<=2&& eta2>-2&& eta2<=-1.5) return 0.968254;
  if(eta1 >1.5&& eta1<=2&& eta2>-1.5&& eta2<=-1) return 0.987764;
  if(eta1 >1.5&& eta1<=2&& eta2>-1&& eta2<=-0.5) return 0.989606;
  if(eta1 >1.5&& eta1<=2&& eta2>-0.5&& eta2<=0) return 0.990423;
  if(eta1 >1.5&& eta1<=2&& eta2>0&& eta2<=0.5) return 0.989735;
  if(eta1 >1.5&& eta1<=2&& eta2>0.5&& eta2<=1) return 0.990882;
  if(eta1 >1.5&& eta1<=2&& eta2>1&& eta2<=1.5) return 0.987764;
  if(eta1 >1.5&& eta1<=2&& eta2>1.5&& eta2<=2) return 0.983886;
  if(eta1 >2&& eta1<=2.5&& eta2>-2.5&& eta2<=-2) return 1.;
  if(eta1 >2&& eta1<=2.5&& eta2>-2&& eta2<=-1.5) return 0.935484;
  if(eta1 >2&& eta1<=2.5&& eta2>-1.5&& eta2<=-1) return 1.;
  if(eta1 >2&& eta1<=2.5&& eta2>-1&& eta2<=-0.5) return 0.984919;
  if(eta1 >2&& eta1<=2.5&& eta2>-0.5&& eta2<=0) return 0.982533;
  if(eta1 >2&& eta1<=2.5&& eta2>0&& eta2<=0.5) return 0.985398;
  if(eta1 >2&& eta1<=2.5&& eta2>0.5&& eta2<=1) return 0.985532;
  if(eta1 >2&& eta1<=2.5&& eta2>1&& eta2<=1.5) return 0.984224;
  if(eta1 >2&& eta1<=2.5&& eta2>1.5&& eta2<=2) return 0.978444;
  if(eta1 >2&& eta1<=2.5&& eta2>2&& eta2<=2.5) return 0.970051;
  return 0.;
}




double Mu17_Mu8orTkMu8_LeadingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 2.58611e-05;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.25) return 0.00104957;
  if(pt >5&& pt<=6&& abseta>1.25&& abseta<=1.6) return 6.7666e-06;
  if(pt >5&& pt<=6&& abseta>1.6&& abseta<=2.1) return 4.84551e-06;
  if(pt >5&& pt<=6&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0.000180462;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.25) return 0;
  if(pt >6&& pt<=7&& abseta>1.25&& abseta<=1.6) return 0;
  if(pt >6&& pt<=7&& abseta>1.6&& abseta<=2.1) return 0.000325903;
  if(pt >6&& pt<=7&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.000308282;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.25) return 0.000447611;
  if(pt >7&& pt<=10&& abseta>1.25&& abseta<=1.6) return 0.000260087;
  if(pt >7&& pt<=10&& abseta>1.6&& abseta<=2.1) return 1.88684e-05;
  if(pt >7&& pt<=10&& abseta>2.1&& abseta<=2.4) return 3.29407e-05;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.000591933;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.25) return 0.00171302;
  if(pt >10&& pt<=15&& abseta>1.25&& abseta<=1.6) return 0.00163596;
  if(pt >10&& pt<=15&& abseta>1.6&& abseta<=2.1) return 0.000343654;
  if(pt >10&& pt<=15&& abseta>2.1&& abseta<=2.4) return 0.0015432;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.646605;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.25) return 0.63285;
  if(pt >15&& pt<=20&& abseta>1.25&& abseta<=1.6) return 0.636385;
  if(pt >15&& pt<=20&& abseta>1.6&& abseta<=2.1) return 0.614301;
  if(pt >15&& pt<=20&& abseta>2.1&& abseta<=2.4) return 0.571553;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.959262;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.25) return 0.958403;
  if(pt >20&& pt<=25&& abseta>1.25&& abseta<=1.6) return 0.971346;
  if(pt >20&& pt<=25&& abseta>1.6&& abseta<=2.1) return 0.950321;
  if(pt >20&& pt<=25&& abseta>2.1&& abseta<=2.4) return 0.910672;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.961259;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.25) return 0.961889;
  if(pt >25&& pt<=30&& abseta>1.25&& abseta<=1.6) return 0.975423;
  if(pt >25&& pt<=30&& abseta>1.6&& abseta<=2.1) return 0.954375;
  if(pt >25&& pt<=30&& abseta>2.1&& abseta<=2.4) return 0.93116;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.960847;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.25) return 0.961149;
  if(pt >30&& pt<=35&& abseta>1.25&& abseta<=1.6) return 0.97578;
  if(pt >30&& pt<=35&& abseta>1.6&& abseta<=2.1) return 0.955396;
  if(pt >30&& pt<=35&& abseta>2.1&& abseta<=2.4) return 0.938805;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.96001;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.25) return 0.960493;
  if(pt >35&& pt<=40&& abseta>1.25&& abseta<=1.6) return 0.975688;
  if(pt >35&& pt<=40&& abseta>1.6&& abseta<=2.1) return 0.955281;
  if(pt >35&& pt<=40&& abseta>2.1&& abseta<=2.4) return 0.943635;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.959531;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.25) return 0.960545;
  if(pt >40&& pt<=50&& abseta>1.25&& abseta<=1.6) return 0.974908;
  if(pt >40&& pt<=50&& abseta>1.6&& abseta<=2.1) return 0.9553;
  if(pt >40&& pt<=50&& abseta>2.1&& abseta<=2.4) return 0.94839;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.958236;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.25) return 0.960457;
  if(pt >50&& pt<=75&& abseta>1.25&& abseta<=1.6) return 0.973904;
  if(pt >50&& pt<=75&& abseta>1.6&& abseta<=2.1) return 0.957003;
  if(pt >50&& pt<=75&& abseta>2.1&& abseta<=2.4) return 0.95065;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.960568;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.25) return 0.962242;
  if(pt >75&& pt<=100&& abseta>1.25&& abseta<=1.6) return 0.972611;
  if(pt >75&& pt<=100&& abseta>1.6&& abseta<=2.1) return 0.959385;
  if(pt >75&& pt<=100&& abseta>2.1&& abseta<=2.4) return 0.958795;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.948178;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.25) return 0.954988;
  if(pt >100&& pt<=200&& abseta>1.25&& abseta<=1.6) return 0.970113;
  if(pt >100&& pt<=200&& abseta>1.6&& abseta<=2.1) return 0.956967;
  if(pt >100&& pt<=200&& abseta>2.1&& abseta<=2.4) return 0.951256;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.932923;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.25) return 0.929771;
  if(pt >200&& pt<=300&& abseta>1.25&& abseta<=1.6) return 0.963679;
  if(pt >200&& pt<=300&& abseta>1.6&& abseta<=2.1) return 0.93643;
  if(pt >200&& pt<=300&& abseta>2.1&& abseta<=2.4) return 0.975752;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.912684;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.25) return 0.940115;
  if(pt >300&& pt<=500&& abseta>1.25&& abseta<=1.6) return 0.960418;
  if(pt >300&& pt<=500&& abseta>1.6&& abseta<=2.1) return 0.955205;
  if(pt >300&& pt<=500&& abseta>2.1&& abseta<=2.4) return 1;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.912684;
  if(pt >500&& abseta>0.8&& abseta<=1.25) return 0.940115;
  if(pt >500&& abseta>1.25&& abseta<=1.6) return 0.960418;
  if(pt >500&& abseta>1.6&& abseta<=2.1) return 0.955205;
  if(pt >500&& abseta>2.1&& abseta<=2.4) return 1;
  return 0;

}




double Mu17_Mu8orTkMu8_TrailingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
if(pt<5) return 0; 
if(abseta<0) return 0; 
if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0.00179503;
if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.25) return 0.0130437;
if(pt >5&& pt<=6&& abseta>1.25&& abseta<=1.6) return 0.00747499;
if(pt >5&& pt<=6&& abseta>1.6&& abseta<=2.1) return 0.00684156;
if(pt >5&& pt<=6&& abseta>2.1&& abseta<=2.4) return 0.0117073;
if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0.00298072;
if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.25) return 0.0181713;
if(pt >6&& pt<=7&& abseta>1.25&& abseta<=1.6) return 0.00965066;
if(pt >6&& pt<=7&& abseta>1.6&& abseta<=2.1) return 0.0100015;
if(pt >6&& pt<=7&& abseta>2.1&& abseta<=2.4) return 0.0167378;
if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.593344;
if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.25) return 0.635037;
if(pt >7&& pt<=10&& abseta>1.25&& abseta<=1.6) return 0.655086;
if(pt >7&& pt<=10&& abseta>1.6&& abseta<=2.1) return 0.631612;
if(pt >7&& pt<=10&& abseta>2.1&& abseta<=2.4) return 0.626282;
if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.934467;
if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.25) return 0.957719;
if(pt >10&& pt<=15&& abseta>1.25&& abseta<=1.6) return 0.982907;
if(pt >10&& pt<=15&& abseta>1.6&& abseta<=2.1) return 0.946048;
if(pt >10&& pt<=15&& abseta>2.1&& abseta<=2.4) return 0.942705;
if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.938632;
if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.25) return 0.960706;
if(pt >15&& pt<=20&& abseta>1.25&& abseta<=1.6) return 0.98549;
if(pt >15&& pt<=20&& abseta>1.6&& abseta<=2.1) return 0.947656;
if(pt >15&& pt<=20&& abseta>2.1&& abseta<=2.4) return 0.946853;
if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.941445;
if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.25) return 0.960256;
if(pt >20&& pt<=25&& abseta>1.25&& abseta<=1.6) return 0.986159;
if(pt >20&& pt<=25&& abseta>1.6&& abseta<=2.1) return 0.948734;
if(pt >20&& pt<=25&& abseta>2.1&& abseta<=2.4) return 0.951101;
if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.945241;
if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.25) return 0.961835;
if(pt >25&& pt<=30&& abseta>1.25&& abseta<=1.6) return 0.986908;
if(pt >25&& pt<=30&& abseta>1.6&& abseta<=2.1) return 0.952527;
if(pt >25&& pt<=30&& abseta>2.1&& abseta<=2.4) return 0.955002;
if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.945182;
if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.25) return 0.960172;
if(pt >30&& pt<=35&& abseta>1.25&& abseta<=1.6) return 0.986761;
if(pt >30&& pt<=35&& abseta>1.6&& abseta<=2.1) return 0.95079;
if(pt >30&& pt<=35&& abseta>2.1&& abseta<=2.4) return 0.955675;
if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.944333;
if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.25) return 0.958681;
if(pt >35&& pt<=40&& abseta>1.25&& abseta<=1.6) return 0.986192;
if(pt >35&& pt<=40&& abseta>1.6&& abseta<=2.1) return 0.949503;
if(pt >35&& pt<=40&& abseta>2.1&& abseta<=2.4) return 0.95582;
if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.94312;
if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.25) return 0.958399;
if(pt >40&& pt<=50&& abseta>1.25&& abseta<=1.6) return 0.985507;
if(pt >40&& pt<=50&& abseta>1.6&& abseta<=2.1) return 0.948039;
if(pt >40&& pt<=50&& abseta>2.1&& abseta<=2.4) return 0.955731;
if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.94264;
if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.25) return 0.958022;
if(pt >50&& pt<=75&& abseta>1.25&& abseta<=1.6) return 0.984603;
if(pt >50&& pt<=75&& abseta>1.6&& abseta<=2.1) return 0.948896;
if(pt >50&& pt<=75&& abseta>2.1&& abseta<=2.4) return 0.954504;
if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.945436;
if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.25) return 0.959033;
if(pt >75&& pt<=100&& abseta>1.25&& abseta<=1.6) return 0.984163;
if(pt >75&& pt<=100&& abseta>1.6&& abseta<=2.1) return 0.952654;
if(pt >75&& pt<=100&& abseta>2.1&& abseta<=2.4) return 0.958314;
if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.936778;
if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.25) return 0.953034;
if(pt >100&& pt<=200&& abseta>1.25&& abseta<=1.6) return 0.982835;
if(pt >100&& pt<=200&& abseta>1.6&& abseta<=2.1) return 0.948853;
if(pt >100&& pt<=200&& abseta>2.1&& abseta<=2.4) return 0.956566;
if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.921908;
if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.25) return 0.945794;
if(pt >200&& pt<=300&& abseta>1.25&& abseta<=1.6) return 0.980339;
if(pt >200&& pt<=300&& abseta>1.6&& abseta<=2.1) return 0.94487;
if(pt >200&& pt<=300&& abseta>2.1&& abseta<=2.4) return 0.940653;
if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.902199;
if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.25) return 0.923168;
if(pt >300&& pt<=500&& abseta>1.25&& abseta<=1.6) return 0.97651;
if(pt >300&& pt<=500&& abseta>1.6&& abseta<=2.1) return 0.954545;
if(pt >300&& pt<=500&& abseta>2.1&& abseta<=2.4) return 1;
if(pt >500&& abseta>0&& abseta<=0.8) return 0.902199;
if(pt >500&& abseta>0.8&& abseta<=1.25) return 0.923168;
if(pt >500&& abseta>1.25&& abseta<=1.6) return 0.97651;
if(pt >500&& abseta>1.6&& abseta<=2.1) return 0.954545;
if(pt >500&& abseta>2.1&& abseta<=2.4) return 1;
return 0;
}




double Mu17_Mu8orTkMu8_DZLeg_Eff(double eta1, double eta2){

  if(eta1 < eta2) {
    double buf =eta1;
    eta1=eta2;
    eta2=buf;
  }
  if(abs(eta1)>2.5) return 0.;
  if(abs(eta2)>2.5) return 0.;
  
  if(eta1 >-2.5&& eta1<=-2&& eta2>-2.5&& eta2<=-2) return 0.972707;
  if(eta1 >-2&& eta1<=-1.5&& eta2>-2.5&& eta2<=-2) return 0.975146;
  if(eta1 >-2&& eta1<=-1.5&& eta2>-2&& eta2<=-1.5) return 0.972411;
  if(eta1 >-1.5&& eta1<=-1&& eta2>-2.5&& eta2<=-2) return 0.973371;
  if(eta1 >-1.5&& eta1<=-1&& eta2>-2&& eta2<=-1.5) return 0.972947;
  if(eta1 >-1.5&& eta1<=-1&& eta2>-1.5&& eta2<=-1) return 0.978426;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-2.5&& eta2<=-2) return 0.97995;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-2&& eta2<=-1.5) return 0.981893;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-1.5&& eta2<=-1) return 0.9832;
  if(eta1 >-1&& eta1<=-0.5&& eta2>-1&& eta2<=-0.5) return 0.987504;
  if(eta1 >-0.5&& eta1<=0&& eta2>-2.5&& eta2<=-2) return 0.983911;
  if(eta1 >-0.5&& eta1<=0&& eta2>-2&& eta2<=-1.5) return 0.983727;
  if(eta1 >-0.5&& eta1<=0&& eta2>-1.5&& eta2<=-1) return 0.983577;
  if(eta1 >-0.5&& eta1<=0&& eta2>-1&& eta2<=-0.5) return 0.988804;
  if(eta1 >-0.5&& eta1<=0&& eta2>-0.5&& eta2<=0) return 0.987586;
  if(eta1 >0&& eta1<=0.5&& eta2>-2.5&& eta2<=-2) return 0.977957;
  if(eta1 >0&& eta1<=0.5&& eta2>-2&& eta2<=-1.5) return 0.976627;
  if(eta1 >0&& eta1<=0.5&& eta2>-1.5&& eta2<=-1) return 0.97841;
  if(eta1 >0&& eta1<=0.5&& eta2>-1&& eta2<=-0.5) return 0.982161;
  if(eta1 >0&& eta1<=0.5&& eta2>-0.5&& eta2<=0) return 0.982768;
  if(eta1 >0&& eta1<=0.5&& eta2>0&& eta2<=0.5) return 0.976456;
  if(eta1 >0.5&& eta1<=1&& eta2>-2.5&& eta2<=-2) return 0.981454;
  if(eta1 >0.5&& eta1<=1&& eta2>-2&& eta2<=-1.5) return 0.974367;
  if(eta1 >0.5&& eta1<=1&& eta2>-1.5&& eta2<=-1) return 0.980125;
  if(eta1 >0.5&& eta1<=1&& eta2>-1&& eta2<=-0.5) return 0.983967;
  if(eta1 >0.5&& eta1<=1&& eta2>-0.5&& eta2<=0) return 0.985229;
  if(eta1 >0.5&& eta1<=1&& eta2>0&& eta2<=0.5) return 0.979161;
  if(eta1 >0.5&& eta1<=1&& eta2>0.5&& eta2<=1) return 0.979705;
  if(eta1 >1&& eta1<=1.5&& eta2>-2.5&& eta2<=-2) return 0.9875;
  if(eta1 >1&& eta1<=1.5&& eta2>-2&& eta2<=-1.5) return 0.983796;
  if(eta1 >1&& eta1<=1.5&& eta2>-1.5&& eta2<=-1) return 0.977081;
  if(eta1 >1&& eta1<=1.5&& eta2>-1&& eta2<=-0.5) return 0.984437;
  if(eta1 >1&& eta1<=1.5&& eta2>-0.5&& eta2<=0) return 0.980074;
  if(eta1 >1&& eta1<=1.5&& eta2>0&& eta2<=0.5) return 0.975956;
  if(eta1 >1&& eta1<=1.5&& eta2>0.5&& eta2<=1) return 0.975936;
  if(eta1 >1&& eta1<=1.5&& eta2>1&& eta2<=1.5) return 0.975914;
  if(eta1 >1.5&& eta1<=2&& eta2>-2.5&& eta2<=-2) return 1.;
  if(eta1 >1.5&& eta1<=2&& eta2>-2&& eta2<=-1.5) return 0.984211;
  if(eta1 >1.5&& eta1<=2&& eta2>-1.5&& eta2<=-1) return 0.984632;
  if(eta1 >1.5&& eta1<=2&& eta2>-1&& eta2<=-0.5) return 0.984092;
  if(eta1 >1.5&& eta1<=2&& eta2>-0.5&& eta2<=0) return 0.980701;
  if(eta1 >1.5&& eta1<=2&& eta2>0&& eta2<=0.5) return 0.976597;
  if(eta1 >1.5&& eta1<=2&& eta2>0.5&& eta2<=1) return 0.978571;
  if(eta1 >1.5&& eta1<=2&& eta2>1&& eta2<=1.5) return 0.975928;
  if(eta1 >1.5&& eta1<=2&& eta2>1.5&& eta2<=2) return 0.979104;
  if(eta1 >2&& eta1<=2.5&& eta2>-2.5&& eta2<=-2) return 0.8;
  if(eta1 >2&& eta1<=2.5&& eta2>-2&& eta2<=-1.5) return 0.979167;
  if(eta1 >2&& eta1<=2.5&& eta2>-1.5&& eta2<=-1) return 0.977823;
  if(eta1 >2&& eta1<=2.5&& eta2>-1&& eta2<=-0.5) return 0.988174;
  if(eta1 >2&& eta1<=2.5&& eta2>-0.5&& eta2<=0) return 0.983944;
  if(eta1 >2&& eta1<=2.5&& eta2>0&& eta2<=0.5) return 0.97293;
  if(eta1 >2&& eta1<=2.5&& eta2>0.5&& eta2<=1) return 0.976471;
  if(eta1 >2&& eta1<=2.5&& eta2>1&& eta2<=1.5) return 0.974463;
  if(eta1 >2&& eta1<=2.5&& eta2>1.5&& eta2<=2) return 0.976945;
  if(eta1 >2&& eta1<=2.5&& eta2>2&& eta2<=2.5) return 0.977016;
  if(eta1 >2.5&& eta2>-2.5&& eta2<=-2) return 0.8;
  if(eta1 >2.5&& eta2>-2&& eta2<=-1.5) return 0.979167;
  if(eta1 >2.5&& eta2>-1.5&& eta2<=-1) return 0.977823;
  if(eta1 >2.5&& eta2>-1&& eta2<=-0.5) return 0.988174;
  if(eta1 >2.5&& eta2>-0.5&& eta2<=0) return 0.983944;
  if(eta1 >2.5&& eta2>0&& eta2<=0.5) return 0.97293;
  if(eta1 >2.5&& eta2>0.5&& eta2<=1) return 0.976471;
  if(eta1 >2.5&& eta2>1&& eta2<=1.5) return 0.974463;
  if(eta1 >2.5&& eta2>1.5&& eta2<=2) return 0.976945;
  if(eta1 >2.5&& eta2>2&& eta2<=2.5) return 0.977016;
  return 0.;
}


double Mu23_Ele8_LeadingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.25) return 0;
  if(pt >5&& pt<=6&& abseta>1.25&& abseta<=1.6) return 0;
  if(pt >5&& pt<=6&& abseta>1.6&& abseta<=2.1) return 0;
  if(pt >5&& pt<=6&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.25) return 0;
  if(pt >6&& pt<=7&& abseta>1.25&& abseta<=1.6) return 0;
  if(pt >6&& pt<=7&& abseta>1.6&& abseta<=2.1) return 0;
  if(pt >6&& pt<=7&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.000271924;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.25) return 0.000271592;
  if(pt >7&& pt<=10&& abseta>1.25&& abseta<=1.6) return 0.000666667;
  if(pt >7&& pt<=10&& abseta>1.6&& abseta<=2.1) return 0.000290782;
  if(pt >7&& pt<=10&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 9.993e-05;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.25) return 0;
  if(pt >10&& pt<=15&& abseta>1.25&& abseta<=1.6) return 0;
  if(pt >10&& pt<=15&& abseta>1.6&& abseta<=2.1) return 0;
  if(pt >10&& pt<=15&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.000177195;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.25) return 0.00287246;
  if(pt >15&& pt<=20&& abseta>1.25&& abseta<=1.6) return 0.00226116;
  if(pt >15&& pt<=20&& abseta>1.6&& abseta<=2.1) return 0.00113058;
  if(pt >15&& pt<=20&& abseta>2.1&& abseta<=2.4) return 0.00341064;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.36652;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.25) return 0.357424;
  if(pt >20&& pt<=25&& abseta>1.25&& abseta<=1.6) return 0.352909;
  if(pt >20&& pt<=25&& abseta>1.6&& abseta<=2.1) return 0.314177;
  if(pt >20&& pt<=25&& abseta>2.1&& abseta<=2.4) return 0.270861;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.924249;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.25) return 0.894449;
  if(pt >25&& pt<=30&& abseta>1.25&& abseta<=1.6) return 0.888675;
  if(pt >25&& pt<=30&& abseta>1.6&& abseta<=2.1) return 0.817686;
  if(pt >25&& pt<=30&& abseta>2.1&& abseta<=2.4) return 0.695873;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.924151;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.25) return 0.902201;
  if(pt >30&& pt<=35&& abseta>1.25&& abseta<=1.6) return 0.899201;
  if(pt >30&& pt<=35&& abseta>1.6&& abseta<=2.1) return 0.820489;
  if(pt >30&& pt<=35&& abseta>2.1&& abseta<=2.4) return 0.748452;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.92661;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.25) return 0.900935;
  if(pt >35&& pt<=40&& abseta>1.25&& abseta<=1.6) return 0.925291;
  if(pt >35&& pt<=40&& abseta>1.6&& abseta<=2.1) return 0.834464;
  if(pt >35&& pt<=40&& abseta>2.1&& abseta<=2.4) return 0.801661;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.922228;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.25) return 0.908025;
  if(pt >40&& pt<=50&& abseta>1.25&& abseta<=1.6) return 0.911715;
  if(pt >40&& pt<=50&& abseta>1.6&& abseta<=2.1) return 0.841576;
  if(pt >40&& pt<=50&& abseta>2.1&& abseta<=2.4) return 0.781871;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.91922;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.25) return 0.906382;
  if(pt >50&& pt<=75&& abseta>1.25&& abseta<=1.6) return 0.902613;
  if(pt >50&& pt<=75&& abseta>1.6&& abseta<=2.1) return 0.845093;
  if(pt >50&& pt<=75&& abseta>2.1&& abseta<=2.4) return 0.821065;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.913009;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.25) return 0.898023;
  if(pt >75&& pt<=100&& abseta>1.25&& abseta<=1.6) return 0.866882;
  if(pt >75&& pt<=100&& abseta>1.6&& abseta<=2.1) return 0.842715;
  if(pt >75&& pt<=100&& abseta>2.1&& abseta<=2.4) return 0.813704;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.912558;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.25) return 0.885558;
  if(pt >100&& pt<=200&& abseta>1.25&& abseta<=1.6) return 0.839119;
  if(pt >100&& pt<=200&& abseta>1.6&& abseta<=2.1) return 0.837785;
  if(pt >100&& pt<=200&& abseta>2.1&& abseta<=2.4) return 0.810109;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.875607;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.25) return 0.854503;
  if(pt >200&& pt<=300&& abseta>1.25&& abseta<=1.6) return 0.811966;
  if(pt >200&& pt<=300&& abseta>1.6&& abseta<=2.1) return 0.842697;
  if(pt >200&& pt<=300&& abseta>2.1&& abseta<=2.4) return 0.645161;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.858696;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.25) return 0.789474;
  if(pt >300&& pt<=500&& abseta>1.25&& abseta<=1.6) return 0.666667;
  if(pt >300&& pt<=500&& abseta>1.6&& abseta<=2.1) return 0.913043;
  if(pt >300&& pt<=500&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.858696;
  if(pt >500&& abseta>0.8&& abseta<=1.25) return 0.789474;
  if(pt >500&& abseta>1.25&& abseta<=1.6) return 0.666667;
  if(pt >500&& abseta>1.6&& abseta<=2.1) return 0.913043;
  if(pt >500&& abseta>2.1&& abseta<=2.4) return 0;
  return 0;


}





double Mu23_Ele8_TrailingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >5&& pt<=6&& abseta>1.479&& abseta<=2) return 0;
  if(pt >5&& pt<=6&& abseta>2&& abseta<=2.5) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >6&& pt<=7&& abseta>1.479&& abseta<=2) return 0;
  if(pt >6&& pt<=7&& abseta>2&& abseta<=2.5) return 0.00362319;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.621256;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.479) return 0.591241;
  if(pt >7&& pt<=10&& abseta>1.479&& abseta<=2) return 0.370487;
  if(pt >7&& pt<=10&& abseta>2&& abseta<=2.5) return 0.222462;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.890974;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.479) return 0.907205;
  if(pt >10&& pt<=15&& abseta>1.479&& abseta<=2) return 0.872526;
  if(pt >10&& pt<=15&& abseta>2&& abseta<=2.5) return 0.735547;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.927847;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.479) return 0.93418;
  if(pt >15&& pt<=20&& abseta>1.479&& abseta<=2) return 0.955895;
  if(pt >15&& pt<=20&& abseta>2&& abseta<=2.5) return 0.848804;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.936082;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.479) return 0.953308;
  if(pt >20&& pt<=25&& abseta>1.479&& abseta<=2) return 0.975988;
  if(pt >20&& pt<=25&& abseta>2&& abseta<=2.5) return 0.926256;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.942757;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.479) return 0.954987;
  if(pt >25&& pt<=30&& abseta>1.479&& abseta<=2) return 0.980809;
  if(pt >25&& pt<=30&& abseta>2&& abseta<=2.5) return 0.95;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.949137;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.479) return 0.963101;
  if(pt >30&& pt<=35&& abseta>1.479&& abseta<=2) return 0.986769;
  if(pt >30&& pt<=35&& abseta>2&& abseta<=2.5) return 0.963224;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.958949;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.479) return 0.966597;
  if(pt >35&& pt<=40&& abseta>1.479&& abseta<=2) return 0.980852;
  if(pt >35&& pt<=40&& abseta>2&& abseta<=2.5) return 0.988196;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.959079;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.479) return 0.968798;
  if(pt >40&& pt<=50&& abseta>1.479&& abseta<=2) return 0.987195;
  if(pt >40&& pt<=50&& abseta>2&& abseta<=2.5) return 0.978207;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.969907;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.479) return 0.972683;
  if(pt >50&& pt<=75&& abseta>1.479&& abseta<=2) return 0.991313;
  if(pt >50&& pt<=75&& abseta>2&& abseta<=2.5) return 0.990315;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.974586;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.479) return 0.977685;
  if(pt >75&& pt<=100&& abseta>1.479&& abseta<=2) return 0.99417;
  if(pt >75&& pt<=100&& abseta>2&& abseta<=2.5) return 0.994382;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.975794;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.479) return 0.982734;
  if(pt >100&& pt<=200&& abseta>1.479&& abseta<=2) return 0.993397;
  if(pt >100&& pt<=200&& abseta>2&& abseta<=2.5) return 0.992071;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.97148;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.479) return 0.976331;
  if(pt >200&& pt<=300&& abseta>1.479&& abseta<=2) return 1;
  if(pt >200&& pt<=300&& abseta>2&& abseta<=2.5) return 0.984127;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.964072;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.479) return 0.973333;
  if(pt >300&& pt<=500&& abseta>1.479&& abseta<=2) return 1;
  if(pt >300&& pt<=500&& abseta>2&& abseta<=2.5) return 1;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.964072;
  if(pt >500&& abseta>0.8&& abseta<=1.479) return 0.973333;
  if(pt >500&& abseta>1.479&& abseta<=2) return 1;
  if(pt >500&& abseta>2&& abseta<=2.5) return 1;
  return 0;
 
}


double Mu8_Ele23_LeadingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >5&& pt<=6&& abseta>1.479&& abseta<=2) return 0;
  if(pt >5&& pt<=6&& abseta>2&& abseta<=2.5) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >6&& pt<=7&& abseta>1.479&& abseta<=2) return 0;
  if(pt >6&& pt<=7&& abseta>2&& abseta<=2.5) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >7&& pt<=10&& abseta>1.479&& abseta<=2) return 0;
  if(pt >7&& pt<=10&& abseta>2&& abseta<=2.5) return 0;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.000136968;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >10&& pt<=15&& abseta>1.479&& abseta<=2) return 0;
  if(pt >10&& pt<=15&& abseta>2&& abseta<=2.5) return 0;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.000669856;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.479) return 0.000154871;
  if(pt >15&& pt<=20&& abseta>1.479&& abseta<=2) return 0.00349345;
  if(pt >15&& pt<=20&& abseta>2&& abseta<=2.5) return 0;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.368521;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.479) return 0.325938;
  if(pt >20&& pt<=25&& abseta>1.479&& abseta<=2) return 0.359808;
  if(pt >20&& pt<=25&& abseta>2&& abseta<=2.5) return 0.230956;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.924471;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.479) return 0.926113;
  if(pt >25&& pt<=30&& abseta>1.479&& abseta<=2) return 0.939342;
  if(pt >25&& pt<=30&& abseta>2&& abseta<=2.5) return 0.879104;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.938433;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.479) return 0.950979;
  if(pt >30&& pt<=35&& abseta>1.479&& abseta<=2) return 0.977019;
  if(pt >30&& pt<=35&& abseta>2&& abseta<=2.5) return 0.93349;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.952497;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.479) return 0.959768;
  if(pt >35&& pt<=40&& abseta>1.479&& abseta<=2) return 0.975772;
  if(pt >35&& pt<=40&& abseta>2&& abseta<=2.5) return 0.973862;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.956532;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.479) return 0.96601;
  if(pt >40&& pt<=50&& abseta>1.479&& abseta<=2) return 0.984054;
  if(pt >40&& pt<=50&& abseta>2&& abseta<=2.5) return 0.968796;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.968399;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.479) return 0.971086;
  if(pt >50&& pt<=75&& abseta>1.479&& abseta<=2) return 0.989417;
  if(pt >50&& pt<=75&& abseta>2&& abseta<=2.5) return 0.984089;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.973664;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.479) return 0.976372;
  if(pt >75&& pt<=100&& abseta>1.479&& abseta<=2) return 0.992112;
  if(pt >75&& pt<=100&& abseta>2&& abseta<=2.5) return 0.988764;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.974848;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.479) return 0.982163;
  if(pt >100&& pt<=200&& abseta>1.479&& abseta<=2) return 0.991746;
  if(pt >100&& pt<=200&& abseta>2&& abseta<=2.5) return 0.99108;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.969697;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.479) return 0.976331;
  if(pt >200&& pt<=300&& abseta>1.479&& abseta<=2) return 1;
  if(pt >200&& pt<=300&& abseta>2&& abseta<=2.5) return 0.984127;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.964072;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.479) return 0.973333;
  if(pt >300&& pt<=500&& abseta>1.479&& abseta<=2) return 1;
  if(pt >300&& pt<=500&& abseta>2&& abseta<=2.5) return 1;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.964072;
  if(pt >500&& abseta>0.8&& abseta<=1.479) return 0.973333;
  if(pt >500&& abseta>1.479&& abseta<=2) return 1;
  if(pt >500&& abseta>2&& abseta<=2.5) return 1;
  return 0;
 
  
}
double Mu8_Ele23_TrailingLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0.000849308;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.25) return 0.00337306;
  if(pt >5&& pt<=6&& abseta>1.25&& abseta<=1.6) return 0.00264996;
  if(pt >5&& pt<=6&& abseta>1.6&& abseta<=2.1) return 0.000556174;
  if(pt >5&& pt<=6&& abseta>2.1&& abseta<=2.4) return 0.0012826;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0.000953516;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.25) return 0.0066762;
  if(pt >6&& pt<=7&& abseta>1.25&& abseta<=1.6) return 0.00319489;
  if(pt >6&& pt<=7&& abseta>1.6&& abseta<=2.1) return 0.00349803;
  if(pt >6&& pt<=7&& abseta>2.1&& abseta<=2.4) return 0.000953289;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.550238;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.25) return 0.56138;
  if(pt >7&& pt<=10&& abseta>1.25&& abseta<=1.6) return 0.551;
  if(pt >7&& pt<=10&& abseta>1.6&& abseta<=2.1) return 0.456819;
  if(pt >7&& pt<=10&& abseta>2.1&& abseta<=2.4) return 0.372804;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.919656;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.25) return 0.929271;
  if(pt >10&& pt<=15&& abseta>1.25&& abseta<=1.6) return 0.943619;
  if(pt >10&& pt<=15&& abseta>1.6&& abseta<=2.1) return 0.837749;
  if(pt >10&& pt<=15&& abseta>2.1&& abseta<=2.4) return 0.783109;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.932312;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.25) return 0.933167;
  if(pt >15&& pt<=20&& abseta>1.25&& abseta<=1.6) return 0.949124;
  if(pt >15&& pt<=20&& abseta>1.6&& abseta<=2.1) return 0.868852;
  if(pt >15&& pt<=20&& abseta>2.1&& abseta<=2.4) return 0.83015;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.931911;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.25) return 0.929875;
  if(pt >20&& pt<=25&& abseta>1.25&& abseta<=1.6) return 0.944263;
  if(pt >20&& pt<=25&& abseta>1.6&& abseta<=2.1) return 0.864733;
  if(pt >20&& pt<=25&& abseta>2.1&& abseta<=2.4) return 0.854305;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.932847;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.25) return 0.928106;
  if(pt >25&& pt<=30&& abseta>1.25&& abseta<=1.6) return 0.94365;
  if(pt >25&& pt<=30&& abseta>1.6&& abseta<=2.1) return 0.871452;
  if(pt >25&& pt<=30&& abseta>2.1&& abseta<=2.4) return 0.834902;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.930261;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.25) return 0.923828;
  if(pt >30&& pt<=35&& abseta>1.25&& abseta<=1.6) return 0.935157;
  if(pt >30&& pt<=35&& abseta>1.6&& abseta<=2.1) return 0.862385;
  if(pt >30&& pt<=35&& abseta>2.1&& abseta<=2.4) return 0.847523;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.931845;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.25) return 0.922326;
  if(pt >35&& pt<=40&& abseta>1.25&& abseta<=1.6) return 0.948252;
  if(pt >35&& pt<=40&& abseta>1.6&& abseta<=2.1) return 0.866011;
  if(pt >35&& pt<=40&& abseta>2.1&& abseta<=2.4) return 0.878229;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.928169;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.25) return 0.929911;
  if(pt >40&& pt<=50&& abseta>1.25&& abseta<=1.6) return 0.9341;
  if(pt >40&& pt<=50&& abseta>1.6&& abseta<=2.1) return 0.868979;
  if(pt >40&& pt<=50&& abseta>2.1&& abseta<=2.4) return 0.843275;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.925363;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.25) return 0.927489;
  if(pt >50&& pt<=75&& abseta>1.25&& abseta<=1.6) return 0.922993;
  if(pt >50&& pt<=75&& abseta>1.6&& abseta<=2.1) return 0.868947;
  if(pt >50&& pt<=75&& abseta>2.1&& abseta<=2.4) return 0.870679;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.918135;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.25) return 0.919188;
  if(pt >75&& pt<=100&& abseta>1.25&& abseta<=1.6) return 0.886571;
  if(pt >75&& pt<=100&& abseta>1.6&& abseta<=2.1) return 0.867803;
  if(pt >75&& pt<=100&& abseta>2.1&& abseta<=2.4) return 0.857602;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.918393;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.25) return 0.910068;
  if(pt >100&& pt<=200&& abseta>1.25&& abseta<=1.6) return 0.862483;
  if(pt >100&& pt<=200&& abseta>1.6&& abseta<=2.1) return 0.869599;
  if(pt >100&& pt<=200&& abseta>2.1&& abseta<=2.4) return 0.862022;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.888241;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.25) return 0.882217;
  if(pt >200&& pt<=300&& abseta>1.25&& abseta<=1.6) return 0.820513;
  if(pt >200&& pt<=300&& abseta>1.6&& abseta<=2.1) return 0.876404;
  if(pt >200&& pt<=300&& abseta>2.1&& abseta<=2.4) return 0.741935;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.86413;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.25) return 0.815789;
  if(pt >300&& pt<=500&& abseta>1.25&& abseta<=1.6) return 0.787879;
  if(pt >300&& pt<=500&& abseta>1.6&& abseta<=2.1) return 0.956522;
  if(pt >300&& pt<=500&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.86413;
  if(pt >500&& abseta>0.8&& abseta<=1.25) return 0.815789;
  if(pt >500&& abseta>1.25&& abseta<=1.6) return 0.787879;
  if(pt >500&& abseta>1.6&& abseta<=2.1) return 0.956522;
  if(pt >500&& abseta>2.1&& abseta<=2.4) return 0;
  return 0;

}

double MuX_EleY_DZLeg_Eff(double eta_e, double eta_mu){
  if(abs(eta_e)>2.5) return 0.;
  if(abs(eta_mu)>2.5) return 0.;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>-2.5&& eta_mu<=-2) return 0.987342;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>-2&& eta_mu<=-1.5) return 0.994118;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>-1.5&& eta_mu<=-1) return 0.994012;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>-1&& eta_mu<=-0.5) return 0.993378;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>-0.5&& eta_mu<=0) return 0.993007;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>0&& eta_mu<=0.5) return 0.988889;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>0.5&& eta_mu<=1) return 0.983333;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>1&& eta_mu<=1.5) return 0.964286;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>1.5&& eta_mu<=2) return 1.;
  if(eta_e >-2.5&& eta_e<=-2&& eta_mu>2&& eta_mu<=2.5) return 1.;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>-2.5&& eta_mu<=-2) return 0.987179;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>-2&& eta_mu<=-1.5) return 0.986885;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>-1.5&& eta_mu<=-1) return 0.992084;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>-1&& eta_mu<=-0.5) return 0.997783;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>-0.5&& eta_mu<=0) return 0.99726;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>0&& eta_mu<=0.5) return 0.996377;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>0.5&& eta_mu<=1) return 0.975;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>1&& eta_mu<=1.5) return 1.;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>1.5&& eta_mu<=2) return 0.979167;
  if(eta_e >-2&& eta_e<=-1.5&& eta_mu>2&& eta_mu<=2.5) return 1.;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>-2.5&& eta_mu<=-2) return 0.992933;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>-2&& eta_mu<=-1.5) return 0.998185;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>-1.5&& eta_mu<=-1) return 0.991549;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>-1&& eta_mu<=-0.5) return 0.997519;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>-0.5&& eta_mu<=0) return 0.996241;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>0&& eta_mu<=0.5) return 1.;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>0.5&& eta_mu<=1) return 0.99789;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>1&& eta_mu<=1.5) return 0.986799;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>1.5&& eta_mu<=2) return 0.993197;
  if(eta_e >-1.5&& eta_e<=-1&& eta_mu>2&& eta_mu<=2.5) return 1.;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>-2&& eta_mu<=-1.5) return 0.998299;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>-1.5&& eta_mu<=-1) return 0.995749;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>-1&& eta_mu<=-0.5) return 0.993833;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>-0.5&& eta_mu<=0) return 0.99654;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>0&& eta_mu<=0.5) return 0.995353;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>0.5&& eta_mu<=1) return 0.997436;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>1&& eta_mu<=1.5) return 0.993068;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>1.5&& eta_mu<=2) return 0.996667;
  if(eta_e >-1&& eta_e<=-0.5&& eta_mu>2&& eta_mu<=2.5) return 0.980392;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>-2&& eta_mu<=-1.5) return 1.;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>-1.5&& eta_mu<=-1) return 0.992891;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>-1&& eta_mu<=-0.5) return 0.998296;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>-0.5&& eta_mu<=0) return 0.998324;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>0&& eta_mu<=0.5) return 0.998315;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>0.5&& eta_mu<=1) return 0.998131;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>1&& eta_mu<=1.5) return 0.994723;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>1.5&& eta_mu<=2) return 1.;
  if(eta_e >-0.5&& eta_e<=0&& eta_mu>2&& eta_mu<=2.5) return 1.;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>-2&& eta_mu<=-1.5) return 1.;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>-1.5&& eta_mu<=-1) return 0.991881;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>-1&& eta_mu<=-0.5) return 0.998049;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>-0.5&& eta_mu<=0) return 0.998329;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>0&& eta_mu<=0.5) return 0.998379;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>0.5&& eta_mu<=1) return 0.997396;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>1&& eta_mu<=1.5) return 0.996591;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>1.5&& eta_mu<=2) return 0.998024;
  if(eta_e >0&& eta_e<=0.5&& eta_mu>2&& eta_mu<=2.5) return 1.;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>-2&& eta_mu<=-1.5) return 1.;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>-1.5&& eta_mu<=-1) return 0.984183;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>-1&& eta_mu<=-0.5) return 0.996355;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>-0.5&& eta_mu<=0) return 0.998016;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>0&& eta_mu<=0.5) return 0.999122;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>0.5&& eta_mu<=1) return 0.998282;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>1&& eta_mu<=1.5) return 0.994919;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>1.5&& eta_mu<=2) return 0.997041;
  if(eta_e >0.5&& eta_e<=1&& eta_mu>2&& eta_mu<=2.5) return 0.996324;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>-2&& eta_mu<=-1.5) return 1.;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>-1.5&& eta_mu<=-1) return 0.996743;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>-1&& eta_mu<=-0.5) return 0.997912;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>-0.5&& eta_mu<=0) return 0.998496;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>0&& eta_mu<=0.5) return 0.996099;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>0.5&& eta_mu<=1) return 0.990885;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>1&& eta_mu<=1.5) return 0.997059;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>1.5&& eta_mu<=2) return 0.99211;
  if(eta_e >1&& eta_e<=1.5&& eta_mu>2&& eta_mu<=2.5) return 0.991701;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>-2&& eta_mu<=-1.5) return 1.;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>-1.5&& eta_mu<=-1) return 0.99115;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>-1&& eta_mu<=-0.5) return 0.993902;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>-0.5&& eta_mu<=0) return 0.992218;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>0&& eta_mu<=0.5) return 0.990909;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>0.5&& eta_mu<=1) return 0.988208;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>1&& eta_mu<=1.5) return 0.994709;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>1.5&& eta_mu<=2) return 0.996377;
  if(eta_e >1.5&& eta_e<=2&& eta_mu>2&& eta_mu<=2.5) return 0.993939;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>-2&& eta_mu<=-1.5) return 1.;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>-1.5&& eta_mu<=-1) return 1.;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>-1&& eta_mu<=-0.5) return 1.;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>-0.5&& eta_mu<=0) return 1.;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>0&& eta_mu<=0.5) return 0.981982;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>0.5&& eta_mu<=1) return 0.994118;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>1&& eta_mu<=1.5) return 0.993902;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>1.5&& eta_mu<=2) return 1.;
  if(eta_e >2&& eta_e<=2.5&& eta_mu>2&& eta_mu<=2.5) return 0.986111;
  if(eta_e >2.5&& eta_mu>-2.5&& eta_mu<=-2) return 1.;
  if(eta_e >2.5&& eta_mu>-2&& eta_mu<=-1.5) return 1.;
  if(eta_e >2.5&& eta_mu>-1.5&& eta_mu<=-1) return 1.;
  if(eta_e >2.5&& eta_mu>-1&& eta_mu<=-0.5) return 1.;
  if(eta_e >2.5&& eta_mu>-0.5&& eta_mu<=0) return 1.;
  if(eta_e >2.5&& eta_mu>0&& eta_mu<=0.5) return 0.981982;
  if(eta_e >2.5&& eta_mu>0.5&& eta_mu<=1) return 0.994118;
  if(eta_e >2.5&& eta_mu>1&& eta_mu<=1.5) return 0.993902;
  if(eta_e >2.5&& eta_mu>1.5&& eta_mu<=2) return 1.;
  if(eta_e >2.5&& eta_mu>2&& eta_mu<=2.5) return 0.986111;
  return 0.; 
}



double DoubleMu8_PFHT300_MuLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0.00106838;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.25) return 0.00416233;
  if(pt >5&& pt<=6&& abseta>1.25&& abseta<=1.6) return 0.00519481;
  if(pt >5&& pt<=6&& abseta>1.6&& abseta<=2.1) return 0.00501505;
  if(pt >5&& pt<=6&& abseta>2.1&& abseta<=2.4) return 0.00223214;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0.000619579;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.25) return 0.00369914;
  if(pt >6&& pt<=7&& abseta>1.25&& abseta<=1.6) return 0;
  if(pt >6&& pt<=7&& abseta>1.6&& abseta<=2.1) return 0.0028169;
  if(pt >6&& pt<=7&& abseta>2.1&& abseta<=2.4) return 0.00330033;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.628881;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.25) return 0.639278;
  if(pt >7&& pt<=10&& abseta>1.25&& abseta<=1.6) return 0.634969;
  if(pt >7&& pt<=10&& abseta>1.6&& abseta<=2.1) return 0.575011;
  if(pt >7&& pt<=10&& abseta>2.1&& abseta<=2.4) return 0.5158;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.920283;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.25) return 0.9202;
  if(pt >10&& pt<=15&& abseta>1.25&& abseta<=1.6) return 0.934516;
  if(pt >10&& pt<=15&& abseta>1.6&& abseta<=2.1) return 0.855689;
  if(pt >10&& pt<=15&& abseta>2.1&& abseta<=2.4) return 0.843246;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.928721;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.25) return 0.922366;
  if(pt >15&& pt<=20&& abseta>1.25&& abseta<=1.6) return 0.946278;
  if(pt >15&& pt<=20&& abseta>1.6&& abseta<=2.1) return 0.869166;
  if(pt >15&& pt<=20&& abseta>2.1&& abseta<=2.4) return 0.850089;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.934192;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.25) return 0.927627;
  if(pt >20&& pt<=25&& abseta>1.25&& abseta<=1.6) return 0.949875;
  if(pt >20&& pt<=25&& abseta>1.6&& abseta<=2.1) return 0.875152;
  if(pt >20&& pt<=25&& abseta>2.1&& abseta<=2.4) return 0.850077;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.9341;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.25) return 0.921638;
  if(pt >25&& pt<=30&& abseta>1.25&& abseta<=1.6) return 0.947709;
  if(pt >25&& pt<=30&& abseta>1.6&& abseta<=2.1) return 0.869196;
  if(pt >25&& pt<=30&& abseta>2.1&& abseta<=2.4) return 0.8553;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.935362;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.25) return 0.922114;
  if(pt >30&& pt<=35&& abseta>1.25&& abseta<=1.6) return 0.945769;
  if(pt >30&& pt<=35&& abseta>1.6&& abseta<=2.1) return 0.876552;
  if(pt >30&& pt<=35&& abseta>2.1&& abseta<=2.4) return 0.86019;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.93497;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.25) return 0.919236;
  if(pt >35&& pt<=40&& abseta>1.25&& abseta<=1.6) return 0.950049;
  if(pt >35&& pt<=40&& abseta>1.6&& abseta<=2.1) return 0.870774;
  if(pt >35&& pt<=40&& abseta>2.1&& abseta<=2.4) return 0.857956;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.932803;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.25) return 0.922988;
  if(pt >40&& pt<=50&& abseta>1.25&& abseta<=1.6) return 0.944914;
  if(pt >40&& pt<=50&& abseta>1.6&& abseta<=2.1) return 0.874174;
  if(pt >40&& pt<=50&& abseta>2.1&& abseta<=2.4) return 0.861485;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.931457;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.25) return 0.919369;
  if(pt >50&& pt<=75&& abseta>1.25&& abseta<=1.6) return 0.942392;
  if(pt >50&& pt<=75&& abseta>1.6&& abseta<=2.1) return 0.874379;
  if(pt >50&& pt<=75&& abseta>2.1&& abseta<=2.4) return 0.862171;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.925765;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.25) return 0.91279;
  if(pt >75&& pt<=100&& abseta>1.25&& abseta<=1.6) return 0.937851;
  if(pt >75&& pt<=100&& abseta>1.6&& abseta<=2.1) return 0.873437;
  if(pt >75&& pt<=100&& abseta>2.1&& abseta<=2.4) return 0.869183;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.92252;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.25) return 0.903452;
  if(pt >100&& pt<=200&& abseta>1.25&& abseta<=1.6) return 0.931099;
  if(pt >100&& pt<=200&& abseta>1.6&& abseta<=2.1) return 0.878222;
  if(pt >100&& pt<=200&& abseta>2.1&& abseta<=2.4) return 0.875661;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.904254;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.25) return 0.880178;
  if(pt >200&& pt<=300&& abseta>1.25&& abseta<=1.6) return 0.930403;
  if(pt >200&& pt<=300&& abseta>1.6&& abseta<=2.1) return 0.852373;
  if(pt >200&& pt<=300&& abseta>2.1&& abseta<=2.4) return 0.779528;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.878072;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.25) return 0.855305;
  if(pt >300&& pt<=500&& abseta>1.25&& abseta<=1.6) return 0.925926;
  if(pt >300&& pt<=500&& abseta>1.6&& abseta<=2.1) return 0.848837;
  if(pt >300&& pt<=500&& abseta>2.1&& abseta<=2.4) return 0.7;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.878072;
  if(pt >500&& abseta>0.8&& abseta<=1.25) return 0.855305;
  if(pt >500&& abseta>1.25&& abseta<=1.6) return 0.925926;
  if(pt >500&& abseta>1.6&& abseta<=2.1) return 0.848837;
  if(pt >500&& abseta>2.1&& abseta<=2.4) return 0.7;
  return 0;
 
  
}

double DoubleEle8_PFHT300_EleLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >5&& pt<=6&& abseta>1.479&& abseta<=2) return 0;
  if(pt >5&& pt<=6&& abseta>2&& abseta<=2.5) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >6&& pt<=7&& abseta>1.479&& abseta<=2) return 0;
  if(pt >6&& pt<=7&& abseta>2&& abseta<=2.5) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.456522;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.479) return 0.4;
  if(pt >7&& pt<=10&& abseta>1.479&& abseta<=2) return 0.4;
  if(pt >7&& pt<=10&& abseta>2&& abseta<=2.5) return 0.307692;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.728395;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.479) return 0.787582;
  if(pt >10&& pt<=15&& abseta>1.479&& abseta<=2) return 0.62963;
  if(pt >10&& pt<=15&& abseta>2&& abseta<=2.5) return 0.703704;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.8475;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.479) return 0.899804;
  if(pt >15&& pt<=20&& abseta>1.479&& abseta<=2) return 0.866667;
  if(pt >15&& pt<=20&& abseta>2&& abseta<=2.5) return 0.909091;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.913732;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.479) return 0.926997;
  if(pt >20&& pt<=25&& abseta>1.479&& abseta<=2) return 0.886957;
  if(pt >20&& pt<=25&& abseta>2&& abseta<=2.5) return 0.895238;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.927039;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.479) return 0.94757;
  if(pt >25&& pt<=30&& abseta>1.479&& abseta<=2) return 0.922807;
  if(pt >25&& pt<=30&& abseta>2&& abseta<=2.5) return 0.956835;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.943736;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.479) return 0.950174;
  if(pt >30&& pt<=35&& abseta>1.479&& abseta<=2) return 0.94362;
  if(pt >30&& pt<=35&& abseta>2&& abseta<=2.5) return 0.971591;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.94;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.479) return 0.964245;
  if(pt >35&& pt<=40&& abseta>1.479&& abseta<=2) return 0.964912;
  if(pt >35&& pt<=40&& abseta>2&& abseta<=2.5) return 1;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.963345;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.479) return 0.955901;
  if(pt >40&& pt<=50&& abseta>1.479&& abseta<=2) return 0.963836;
  if(pt >40&& pt<=50&& abseta>2&& abseta<=2.5) return 0.969072;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.957167;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.479) return 0.969603;
  if(pt >50&& pt<=75&& abseta>1.479&& abseta<=2) return 0.980392;
  if(pt >50&& pt<=75&& abseta>2&& abseta<=2.5) return 0.977055;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.962898;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.479) return 0.962046;
  if(pt >75&& pt<=100&& abseta>1.479&& abseta<=2) return 0.969559;
  if(pt >75&& pt<=100&& abseta>2&& abseta<=2.5) return 0.996094;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.951738;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.479) return 0.952117;
  if(pt >100&& pt<=200&& abseta>1.479&& abseta<=2) return 0.966834;
  if(pt >100&& pt<=200&& abseta>2&& abseta<=2.5) return 0.98103;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.930233;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.479) return 0.934459;
  if(pt >200&& pt<=300&& abseta>1.479&& abseta<=2) return 0.913793;
  if(pt >200&& pt<=300&& abseta>2&& abseta<=2.5) return 0.931507;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.950445;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.479) return 0.935574;
  if(pt >300&& pt<=500&& abseta>1.479&& abseta<=2) return 0.935185;
  if(pt >300&& pt<=500&& abseta>2&& abseta<=2.5) return 1;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.950445;
  if(pt >500&& abseta>0.8&& abseta<=1.479) return 0.935574;
  if(pt >500&& abseta>1.479&& abseta<=2) return 0.935185;
  if(pt >500&& abseta>2&& abseta<=2.5) return 1;
  return 0;
  
}

double Mu8Ele8_PFHT300_MuLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.25) return 0;
  if(pt >5&& pt<=6&& abseta>1.25&& abseta<=1.6) return 0;
  if(pt >5&& pt<=6&& abseta>1.6&& abseta<=2.1) return 0;
  if(pt >5&& pt<=6&& abseta>2.1&& abseta<=2.4) return 0.0222222;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.25) return 0;
  if(pt >6&& pt<=7&& abseta>1.25&& abseta<=1.6) return 0;
  if(pt >6&& pt<=7&& abseta>1.6&& abseta<=2.1) return 0;
  if(pt >6&& pt<=7&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.658986;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.25) return 0.593407;
  if(pt >7&& pt<=10&& abseta>1.25&& abseta<=1.6) return 0.655556;
  if(pt >7&& pt<=10&& abseta>1.6&& abseta<=2.1) return 0.591398;
  if(pt >7&& pt<=10&& abseta>2.1&& abseta<=2.4) return 0.568627;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.928425;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.25) return 0.941176;
  if(pt >10&& pt<=15&& abseta>1.25&& abseta<=1.6) return 0.930556;
  if(pt >10&& pt<=15&& abseta>1.6&& abseta<=2.1) return 0.814103;
  if(pt >10&& pt<=15&& abseta>2.1&& abseta<=2.4) return 0.819672;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.92635;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.25) return 0.954545;
  if(pt >15&& pt<=20&& abseta>1.25&& abseta<=1.6) return 0.953368;
  if(pt >15&& pt<=20&& abseta>1.6&& abseta<=2.1) return 0.870466;
  if(pt >15&& pt<=20&& abseta>2.1&& abseta<=2.4) return 0.816901;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.931818;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.25) return 0.936364;
  if(pt >20&& pt<=25&& abseta>1.25&& abseta<=1.6) return 0.935185;
  if(pt >20&& pt<=25&& abseta>1.6&& abseta<=2.1) return 0.853107;
  if(pt >20&& pt<=25&& abseta>2.1&& abseta<=2.4) return 0.851852;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.919622;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.25) return 0.915167;
  if(pt >25&& pt<=30&& abseta>1.25&& abseta<=1.6) return 0.932432;
  if(pt >25&& pt<=30&& abseta>1.6&& abseta<=2.1) return 0.831111;
  if(pt >25&& pt<=30&& abseta>2.1&& abseta<=2.4) return 0.805195;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.925519;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.25) return 0.904908;
  if(pt >30&& pt<=35&& abseta>1.25&& abseta<=1.6) return 0.878641;
  if(pt >30&& pt<=35&& abseta>1.6&& abseta<=2.1) return 0.876238;
  if(pt >30&& pt<=35&& abseta>2.1&& abseta<=2.4) return 0.926471;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.927536;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.25) return 0.904615;
  if(pt >35&& pt<=40&& abseta>1.25&& abseta<=1.6) return 0.917874;
  if(pt >35&& pt<=40&& abseta>1.6&& abseta<=2.1) return 0.890052;
  if(pt >35&& pt<=40&& abseta>2.1&& abseta<=2.4) return 0.848485;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.915;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.25) return 0.931419;
  if(pt >40&& pt<=50&& abseta>1.25&& abseta<=1.6) return 0.931624;
  if(pt >40&& pt<=50&& abseta>1.6&& abseta<=2.1) return 0.879725;
  if(pt >40&& pt<=50&& abseta>2.1&& abseta<=2.4) return 0.87037;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.926584;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.25) return 0.926678;
  if(pt >50&& pt<=75&& abseta>1.25&& abseta<=1.6) return 0.923077;
  if(pt >50&& pt<=75&& abseta>1.6&& abseta<=2.1) return 0.857977;
  if(pt >50&& pt<=75&& abseta>2.1&& abseta<=2.4) return 0.847953;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.907633;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.25) return 0.890601;
  if(pt >75&& pt<=100&& abseta>1.25&& abseta<=1.6) return 0.846354;
  if(pt >75&& pt<=100&& abseta>1.6&& abseta<=2.1) return 0.845361;
  if(pt >75&& pt<=100&& abseta>2.1&& abseta<=2.4) return 0.855422;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.899682;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.25) return 0.886555;
  if(pt >100&& pt<=200&& abseta>1.25&& abseta<=1.6) return 0.840611;
  if(pt >100&& pt<=200&& abseta>1.6&& abseta<=2.1) return 0.8375;
  if(pt >100&& pt<=200&& abseta>2.1&& abseta<=2.4) return 0.780822;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.878613;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.25) return 0.837037;
  if(pt >200&& pt<=300&& abseta>1.25&& abseta<=1.6) return 0.802469;
  if(pt >200&& pt<=300&& abseta>1.6&& abseta<=2.1) return 0.90566;
  if(pt >200&& pt<=300&& abseta>2.1&& abseta<=2.4) return 0.8;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.8;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.25) return 0.744681;
  if(pt >300&& pt<=500&& abseta>1.25&& abseta<=1.6) return 0.722222;
  if(pt >300&& pt<=500&& abseta>1.6&& abseta<=2.1) return 1;
  if(pt >300&& pt<=500&& abseta>2.1&& abseta<=2.4) return 0;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.8;
  if(pt >500&& abseta>0.8&& abseta<=1.25) return 0.744681;
  if(pt >500&& abseta>1.25&& abseta<=1.6) return 0.722222;
  if(pt >500&& abseta>1.6&& abseta<=2.1) return 1;
  if(pt >500&& abseta>2.1&& abseta<=2.4) return 0;
  return 0;
 
  
}

double Mu8Ele8_PFHT300_EleLeg_Eff(double pt, double eta){
  double abseta =abs(eta);
  if(pt<5) return 0; 
  if(abseta<0) return 0; 
  if(pt >5&& pt<=6&& abseta>0&& abseta<=0.8) return 0;
  if(pt >5&& pt<=6&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >5&& pt<=6&& abseta>1.479&& abseta<=2) return 0;
  if(pt >5&& pt<=6&& abseta>2&& abseta<=2.5) return 0;
  if(pt >6&& pt<=7&& abseta>0&& abseta<=0.8) return 0;
  if(pt >6&& pt<=7&& abseta>0.8&& abseta<=1.479) return 0;
  if(pt >6&& pt<=7&& abseta>1.479&& abseta<=2) return 0;
  if(pt >6&& pt<=7&& abseta>2&& abseta<=2.5) return 0;
  if(pt >7&& pt<=10&& abseta>0&& abseta<=0.8) return 0.453947;
  if(pt >7&& pt<=10&& abseta>0.8&& abseta<=1.479) return 0.465;
  if(pt >7&& pt<=10&& abseta>1.479&& abseta<=2) return 0.358974;
  if(pt >7&& pt<=10&& abseta>2&& abseta<=2.5) return 0.206349;
  if(pt >10&& pt<=15&& abseta>0&& abseta<=0.8) return 0.756145;
  if(pt >10&& pt<=15&& abseta>0.8&& abseta<=1.479) return 0.787393;
  if(pt >10&& pt<=15&& abseta>1.479&& abseta<=2) return 0.682759;
  if(pt >10&& pt<=15&& abseta>2&& abseta<=2.5) return 0.75;
  if(pt >15&& pt<=20&& abseta>0&& abseta<=0.8) return 0.882503;
  if(pt >15&& pt<=20&& abseta>0.8&& abseta<=1.479) return 0.901201;
  if(pt >15&& pt<=20&& abseta>1.479&& abseta<=2) return 0.833871;
  if(pt >15&& pt<=20&& abseta>2&& abseta<=2.5) return 0.812785;
  if(pt >20&& pt<=25&& abseta>0&& abseta<=0.8) return 0.923245;
  if(pt >20&& pt<=25&& abseta>0.8&& abseta<=1.479) return 0.950128;
  if(pt >20&& pt<=25&& abseta>1.479&& abseta<=2) return 0.892534;
  if(pt >20&& pt<=25&& abseta>2&& abseta<=2.5) return 0.93038;
  if(pt >25&& pt<=30&& abseta>0&& abseta<=0.8) return 0.935046;
  if(pt >25&& pt<=30&& abseta>0.8&& abseta<=1.479) return 0.954974;
  if(pt >25&& pt<=30&& abseta>1.479&& abseta<=2) return 0.922529;
  if(pt >25&& pt<=30&& abseta>2&& abseta<=2.5) return 0.947126;
  if(pt >30&& pt<=35&& abseta>0&& abseta<=0.8) return 0.953473;
  if(pt >30&& pt<=35&& abseta>0.8&& abseta<=1.479) return 0.95775;
  if(pt >30&& pt<=35&& abseta>1.479&& abseta<=2) return 0.964401;
  if(pt >30&& pt<=35&& abseta>2&& abseta<=2.5) return 0.952965;
  if(pt >35&& pt<=40&& abseta>0&& abseta<=0.8) return 0.961719;
  if(pt >35&& pt<=40&& abseta>0.8&& abseta<=1.479) return 0.964186;
  if(pt >35&& pt<=40&& abseta>1.479&& abseta<=2) return 0.959383;
  if(pt >35&& pt<=40&& abseta>2&& abseta<=2.5) return 0.981405;
  if(pt >40&& pt<=50&& abseta>0&& abseta<=0.8) return 0.964031;
  if(pt >40&& pt<=50&& abseta>0.8&& abseta<=1.479) return 0.971772;
  if(pt >40&& pt<=50&& abseta>1.479&& abseta<=2) return 0.978242;
  if(pt >40&& pt<=50&& abseta>2&& abseta<=2.5) return 0.973176;
  if(pt >50&& pt<=75&& abseta>0&& abseta<=0.8) return 0.971595;
  if(pt >50&& pt<=75&& abseta>0.8&& abseta<=1.479) return 0.972588;
  if(pt >50&& pt<=75&& abseta>1.479&& abseta<=2) return 0.985959;
  if(pt >50&& pt<=75&& abseta>2&& abseta<=2.5) return 0.988869;
  if(pt >75&& pt<=100&& abseta>0&& abseta<=0.8) return 0.977452;
  if(pt >75&& pt<=100&& abseta>0.8&& abseta<=1.479) return 0.978663;
  if(pt >75&& pt<=100&& abseta>1.479&& abseta<=2) return 0.991364;
  if(pt >75&& pt<=100&& abseta>2&& abseta<=2.5) return 0.994388;
  if(pt >100&& pt<=200&& abseta>0&& abseta<=0.8) return 0.980729;
  if(pt >100&& pt<=200&& abseta>0.8&& abseta<=1.479) return 0.984633;
  if(pt >100&& pt<=200&& abseta>1.479&& abseta<=2) return 0.991087;
  if(pt >100&& pt<=200&& abseta>2&& abseta<=2.5) return 0.993088;
  if(pt >200&& pt<=300&& abseta>0&& abseta<=0.8) return 0.978744;
  if(pt >200&& pt<=300&& abseta>0.8&& abseta<=1.479) return 0.981289;
  if(pt >200&& pt<=300&& abseta>1.479&& abseta<=2) return 1;
  if(pt >200&& pt<=300&& abseta>2&& abseta<=2.5) return 1;
  if(pt >300&& pt<=500&& abseta>0&& abseta<=0.8) return 0.968944;
  if(pt >300&& pt<=500&& abseta>0.8&& abseta<=1.479) return 1;
  if(pt >300&& pt<=500&& abseta>1.479&& abseta<=2) return 1;
  if(pt >300&& pt<=500&& abseta>2&& abseta<=2.5) return 1;
  if(pt >500&& abseta>0&& abseta<=0.8) return 0.968944;
  if(pt >500&& abseta>0.8&& abseta<=1.479) return 1;
  if(pt >500&& abseta>1.479&& abseta<=2) return 1;
  if(pt >500&& abseta>2&& abseta<=2.5) return 1;
  return 0;

  
}


double HTLeg_Eff(int lep1_pdgId, double lep1_pt, int lep2_pdgId, double lep2_pt,  double HT){
  if(abs(lep1_pdgId)!=13 || abs(lep2_pdgId) !=13) return 1 ;// HT ineffcy  only measured for dimuon+HT. 
  if(HT<300) return 1.;
  if(HT>1000) return 1.;
  if(HT<1000&&HT>350) return 0.992;//Effect due to the L1 loss of efficiency in Run H  
  if(HT<350 && lep1_pt<27&&lep2_pt<27) return 0.990; //HLT SF in that region  is below 1.
  return 1.;

  
}





double TotalTriggerEff(int lep1_pdgId, double lep1_pt, double lep1_eta, int lep2_pdgId, double lep2_pt, double lep2_eta, double HT){

  double fracrunH =  1-27.22/36.46;
  
  // Make sure lep 1 is the leading lepton
  if(lep1_pt<lep2_pt){
    int buf_pdgId = lep1_pdgId;
    double buf_pt = lep1_pt;
    double buf_eta = lep1_eta;
    lep1_pdgId =  lep2_pdgId;
    lep1_pt = lep2_pt;
    lep1_eta = lep2_eta;
    lep2_pdgId =  buf_pdgId;
    lep2_pt = buf_pt;
    lep2_eta = buf_eta;
    
  }



  if(abs(lep1_pdgId)+ abs(lep2_pdgId) ==22) {//ele-Ele
    if(HT<300) return Ele23_Ele12_LeadingLeg_Eff(lep1_pt,lep1_eta)*Ele23_Ele12_TrailingLeg_Eff(lep2_pt,lep2_eta)*Ele23_Ele12_DZLeg_Eff(lep1_eta,lep2_eta);
    else  return DoubleEle8_PFHT300_EleLeg_Eff(lep1_pt,lep1_eta)*DoubleEle8_PFHT300_EleLeg_Eff(lep2_pt,lep2_eta) ; //HT SF is one here. 
  }
  else if(abs(lep1_pdgId)+ abs(lep2_pdgId) ==26) {//mu-mu
    if(HT<300) return  Mu17_Mu8orTkMu8_LeadingLeg_Eff(lep1_pt,lep1_eta)*Mu17_Mu8orTkMu8_TrailingLeg_Eff(lep2_pt,lep2_eta)*(fracrunH*Mu17_Mu8orTkMu8_DZLeg_Eff(lep1_eta,lep2_eta)+(1-fracrunH));
    //The last factor applies the DZ effcy on the fraction of lumi corresponding to Run H
    else return DoubleMu8_PFHT300_MuLeg_Eff(lep1_pt,lep1_eta)*DoubleMu8_PFHT300_MuLeg_Eff(lep2_pt,lep2_eta)*HTLeg_Eff(lep1_pdgId,lep1_pt,lep2_pdgId,lep2_pt,HT);
  }
  else if(abs(lep1_pdgId)+ abs(lep2_pdgId) ==24) {//ele-mu
    if(abs(lep1_pdgId)==11){//leading ele
      if(HT<300) return Mu8_Ele23_LeadingLeg_Eff(lep1_pt,lep1_eta)*Mu8_Ele23_TrailingLeg_Eff(lep2_pt,lep2_eta)*(fracrunH*MuX_EleY_DZLeg_Eff(lep1_eta,lep2_eta)+(1-fracrunH)); 
      else return Mu8Ele8_PFHT300_EleLeg_Eff(lep1_pt,lep1_eta)*Mu8Ele8_PFHT300_MuLeg_Eff(lep2_pt,lep2_eta);
    }
    else{//leading mu
      if(HT<300) return Mu23_Ele8_LeadingLeg_Eff(lep1_pt,lep1_eta)*Mu23_Ele8_TrailingLeg_Eff(lep2_pt,lep2_eta)*(fracrunH*MuX_EleY_DZLeg_Eff(lep2_eta,lep1_eta)+(1-fracrunH));
      else return Mu8Ele8_PFHT300_EleLeg_Eff(lep2_pt,lep2_eta)*Mu8Ele8_PFHT300_MuLeg_Eff(lep1_pt,lep1_eta);
    } 
  }
  
  else return 0.; 
  
}

double FastSimTriggerEfficiency(int pdgId1, int pdgId2, float pt1, float pt2, float eta1, float eta2, float ht) {
    return TotalTriggerEff(pdgId1, pt1, eta1, pdgId2, pt2, eta2, ht);
}
