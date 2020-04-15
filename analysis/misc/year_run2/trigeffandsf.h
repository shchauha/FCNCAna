#include <map>

bool debug = false;

TFile *myf; // NJA
TString rfilepath = Form("%s/analysis/misc/year_run2/triggeffcymapsRA5_Run2_ALL.root",getenv("FTBASE"));
std::map<TString,TH2D*> hists;

double LegEffcyorSF(TString legname, double pt, double eta, TString year, bool issf, int systfluc){
    if(legname =="ele23_diele" && (year=="2018"||  year=="2017BtoF"||  year=="2017CtoF") ) legname = "ele23l1eg24_diele";
    if(debug)  cout << legname<<"  " ; 
    TString hname = year+"/"+legname;
    TString hname_mc = hname+"_mc";

    TH2D* myhdata; // NJA
    TH2D* myhmc;
    if (!hists.count(hname)) {
        myhdata = (TH2D *)(myf->Get(hname));
        myhmc = (TH2D *)(myf->Get(hname_mc));
        hists[hname] = myhdata;
        hists[hname_mc] = myhmc;
    } else {
        myhdata = hists[hname];
        myhmc = hists[hname_mc];
    }


    int thebin= myhdata->FindBin(pt,eta);
    double effdata = myhdata->GetBinContent(thebin);
    double effmc = myhmc->GetBinContent(thebin);
    if(debug)  cout <<effdata<<endl ; 

    if(systfluc==0&& !issf) return effdata;
    else if(systfluc==0){
        if(effmc >0) return effdata/effmc;
        else return 0;
    }
    else{
        double err_stat = myhdata->GetBinError(thebin);
        double err_syst = 0.01;
        //Special case for dielectron leading leg in 2017/2018 
        if(legname =="ele23l1eg24_diele" && (year=="2018" ||  year=="2017BtoF"||  year=="2017CtoF") ){
            TString alternhistname = (systfluc==1)?  year+"/"+"ele23l1eg22_diele" : year+"/"+"ele23l1eg25_diele"; 

            TH2D* myalthisto; // NJA
            if (!hists.count(alternhistname)) {
                myalthisto =   (TH2D *) myf->Get(alternhistname) ;
                hists[alternhistname] = myalthisto;
            } else {
                myalthisto = hists[alternhistname];
            }

            double var = fabs(myalthisto->GetBinContent(thebin) - effdata); 
            err_syst  = sqrt( pow(err_syst ,2.)+ pow(var,2.));

        }

        double err_tot = sqrt( pow(err_stat,2) + pow(err_syst,2));
        if(systfluc==-1&&!issf ) return TMath::Max(effdata-err_tot,0.);
        if(systfluc==-1&&issf ) return effdata/effmc-err_tot;
        if(systfluc==1&&!issf ) return TMath::Min(effdata+err_tot,1.);
        if(systfluc==1&&issf ) return  effdata/effmc+err_tot ;
    }

    return 1;

}

double HTLeg(bool islowlow, double ht, TString year, bool issf, int systfluc){
    if(year=="2016" ||  year=="2017BtoF"||  year=="2017CtoF") return 1;
    else if(year=="2018" && !islowlow) return 1;
    else if(year=="2018" && islowlow&&ht>450)  return 1.;
    else if(year=="2018" && islowlow&&ht<=450) {
        if(systfluc==0) return 0.95;
        if(systfluc==1) return 1;
        if(systfluc==-1) return 0.9;
    }

    return 1;
}

double DZLeg(bool islowlow, int nmus, TString year, bool issf, int systfluc){

    double effdz_2016[3] ={0.992,0.999,0.995};
    double effdz_2017[3] ={0.991,0.995,0.995};
    double effdz_2018[3] ={0.995,0.997,0.999};

    if(year=="2016" && islowlow) { // No DZ in dilepton+HT triggers in 2016
        return 1.;
    }
    else if( ( year=="2017BtoF"||  year=="2017CtoF" )&& !islowlow && nmus==0) { // No DZ in isolated dielectron triggers in 2017
        return 1.;
    }
    else if(year=="2018" && !islowlow && nmus==0) { // No DZ in isolated dielectron triggers in 2018
        return 1.;
    }
    else if(year=="2016") {
        int fluc = 0;  
        double centralval = (issf)?1. : effdz_2016[nmus]  ;
        if(systfluc==-1) fluc = -1;
        if(systfluc==1) fluc = +1;
        centralval = TMath::Max(TMath::Min(centralval +fluc*0.003, 1.) , 0.);  
        return centralval; 
    }
    else if( year=="2017BtoF"||  year=="2017CtoF" ) {
        int fluc = 0;  
        double centralval = (issf)?1. : effdz_2017[nmus]  ;
        if(systfluc==-1) fluc = -1;
        if(systfluc==1) fluc = +1;
        centralval = TMath::Max(TMath::Min(centralval +fluc*0.003, 1.) , 0.);  
        return centralval; 
    }
    else if(year=="2018") {
        int fluc = 0;  
        double centralval = (issf)?1. : effdz_2018[nmus]  ;
        if(systfluc==-1) fluc = -1;
        if(systfluc==1) fluc = +1;
        centralval = TMath::Max(TMath::Min(centralval +fluc*0.003, 1.) , 0.);  
        return centralval; 
    }

    return 1.;

}

double TriggerWeight(int pdgid1, double pt1, double eta1, int pdgid2, double pt2, double eta2, double ht, int iyear, bool issf, int systfluc){

    if (!myf) { // NJA
        myf = new TFile(rfilepath,"read");
    }

    bool islowlow = pt1<25 && pt2<25 ;

    int nmus (0);
    if(fabs(pdgid1)==13) nmus++;
    if(fabs(pdgid2)==13) nmus++;

    TString year = "2016";
    if (iyear == 2017) year = "2017";
    if (iyear == 2018) year = "2018";


    //Find leading/trailing lepton 
    double ptlead, etalead, pttrail, etatrail;
    int pdgidlead;
    // int pdgidtrail;

    if(pt1>pt2) {
        ptlead = TMath::Min(pt1,499.); 
        etalead = fabs(eta1) ; 
        pdgidlead = pdgid1; 
        pttrail = TMath::Min(pt2,499.);  
        etatrail =fabs(eta2) ; 
        // pdgidtrail = pdgid2; 
    }
    else {
        ptlead =  TMath::Min(pt2,499.);  
        etalead =fabs(eta2) ; 
        pdgidlead = pdgid2; 
        pttrail =TMath::Min(pt1,499.); 
        etatrail =fabs(eta1) ; 
        // pdgidtrail = pdgid1; 
    }


    double leptonlegs=1; 
    double htleg= HTLeg(islowlow, ht, year, issf, systfluc); 
    double dzleg= DZLeg(islowlow, nmus, year, issf, systfluc);
    if(!islowlow &&year=="2017") year = "2017BtoF";
    if(islowlow) {
        if(year=="2017") year = "2017CtoF";
        if(nmus==0)leptonlegs = LegEffcyorSF("ele_dieleht",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("ele_dieleht",pttrail,etatrail,year,issf,systfluc); 
        else if(nmus==1 && fabs(pdgidlead)==13)leptonlegs = LegEffcyorSF("mu_elemuht",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("ele_elemuht",pttrail,etatrail,year,issf,systfluc); 
        else if(nmus==1 && fabs(pdgidlead)==11)leptonlegs = LegEffcyorSF("ele_elemuht",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("mu_elemuht",pttrail,etatrail,year,issf,systfluc); 
        else if(nmus==2)leptonlegs = LegEffcyorSF("mu_elemuht",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("mu_elemuht",pttrail,etatrail,year,issf,systfluc); 
    }
    else{

        if(nmus==0)leptonlegs = LegEffcyorSF("ele23_diele",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("ele12_diele",pttrail,etatrail,year,issf,systfluc); 
        else if(nmus==1 && fabs(pdgidlead)==13)leptonlegs = LegEffcyorSF("mu23_elemu",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("ele12_elemu",pttrail,etatrail,year,issf,systfluc); 
        // below line added to fix bug with 2017BtoF map
        else if(nmus==1 && fabs(pdgidlead)==11&& year == "2017BtoF"){year = "2017CtoF"; leptonlegs = LegEffcyorSF("ele23_elemu",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("mu8_elemu",pttrail,etatrail,year,issf,systfluc); }
        else if(nmus==1 && fabs(pdgidlead)==11)leptonlegs = LegEffcyorSF("ele23_elemu",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("mu8_elemu",pttrail,etatrail,year,issf,systfluc); 
        else if(nmus==2)leptonlegs = LegEffcyorSF("mu17_dimu",ptlead,etalead,year,issf,systfluc) * LegEffcyorSF("mu8_dimu",pttrail,etatrail,year,issf,systfluc); 
    }


    if(debug) cout << "HT leg, dzleg, leptonlegs "<< htleg<<", "<< dzleg <<", "<<leptonlegs <<endl;
    return  htleg*dzleg*leptonlegs;

}




