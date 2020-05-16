#pragma GCC diagnostic ignored "-Wsign-compare"

#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1.h"
#include "TChain.h"

#include <math.h>

#include <signal.h>
#include <stdlib.h>
#include <stdio.h>

#include "CORE/CMS3.h"
#include "CORE/JetSelections.h"
#include "CORE/Config.h"

#include "tqdm.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;
using namespace tas;

bool STOP_REQUESTED = false;

int ScanChain(TChain *ch, TString options, TString outname) {

    signal(SIGINT, [](int){
            cout << "SIGINT Caught, stopping after current event" << endl;
            STOP_REQUESTED=true;
            });
    STOP_REQUESTED=false;

    float WP_DEEPCSV_TIGHT  = -1;
    float WP_DEEPCSV_MEDIUM = -1;
    float WP_DEEPCSV_LOOSE  = -1;
    std::string jecEraMC = "Fall17_17Nov2017_V32";
    if (options.Contains("2016")) {
        WP_DEEPCSV_TIGHT  = 0.8958; // 2016 80X
        WP_DEEPCSV_MEDIUM = 0.6324;
        WP_DEEPCSV_LOOSE  = 0.2219;
        jecEraMC = "Summer16_23Sep2016V4";
    }
    if (options.Contains("2017")) {
        WP_DEEPCSV_TIGHT  = 0.8001; // 2017 94X
        WP_DEEPCSV_MEDIUM = 0.4941;
        WP_DEEPCSV_LOOSE  = 0.1522;
        jecEraMC = "Fall17_17Nov2017_V32";
    }
    if (options.Contains("2018")) {
        WP_DEEPCSV_TIGHT  = 0.7527; // 2018 102X
        WP_DEEPCSV_MEDIUM = 0.4184;
        WP_DEEPCSV_LOOSE  = 0.1241;
        jecEraMC = "Fall17_17Nov2017_V32";
    }
    bool is80X = options.Contains("80X");
    bool quiet = options.Contains("quiet");

    double eps = 0.00001;
    const int Pt_nBins = 17;
    double Pt_min      = 20.0;
    double Pt_max      = 600.0;
    double *Pt_bins    = new double[Pt_nBins+1]{20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 120.0, 150.0, 200.0, 300.0, 400.0, 600.0+eps};
    const int Eta_nBins = 7;
    double Eta_min      = 0.0;
    double Eta_max      = 2.8;
    double *Eta_bins    = new double[Eta_nBins+1];
    for (int iBin=0; iBin<Eta_nBins+1; iBin++){
        Eta_bins[iBin] = (iBin)*(Eta_max-Eta_min)/Eta_nBins;
        if(iBin==Eta_nBins) Eta_bins[iBin] += eps;
    }

    TH2D* h2_BTaggingEff_csv_Denom_b = new TH2D("h2_BTaggingEff_csv_Denom_b", "b-tagging effiency, Denominator, MC, CSV, b-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_Denom_c = new TH2D("h2_BTaggingEff_csv_Denom_c", "b-tagging effiency, Denominator, MC, CSV loose, c-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_Denom_udsg = new TH2D("h2_BTaggingEff_csv_Denom_udsg", "b-tagging effiency, Denominator, MC, CSV loose, udsg-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);

    TH2D* h2_BTaggingEff_csv_tight_Num_b = new TH2D("h2_BTaggingEff_csv_tight_Num_b", "b-tagging effiency, Numerator, MC, CSV tight, b-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_tight_Num_c = new TH2D("h2_BTaggingEff_csv_tight_Num_c", "b-tagging effiency, Numerator, MC, CSV tight, c-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_tight_Num_udsg = new TH2D("h2_BTaggingEff_csv_tight_Num_udsg", "b-tagging effiency, Numerator, MC, CSV tight, udsg-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);

    TH2D* h2_BTaggingEff_csv_med_Num_b = new TH2D("h2_BTaggingEff_csv_med_Num_b", "b-tagging effiency, Numerator, MC, CSV med, b-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_med_Num_c = new TH2D("h2_BTaggingEff_csv_med_Num_c", "b-tagging effiency, Numerator, MC, CSV med, c-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_med_Num_udsg = new TH2D("h2_BTaggingEff_csv_med_Num_udsg", "b-tagging effiency, Numerator, MC, CSV med, udsg-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);

    TH2D* h2_BTaggingEff_csv_loose_Num_b = new TH2D("h2_BTaggingEff_csv_loose_Num_b", "b-tagging effiency, Numerator, MC, CSV loose, b-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_loose_Num_c = new TH2D("h2_BTaggingEff_csv_loose_Num_c", "b-tagging effiency, Numerator, MC, CSV loose, c-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);
    TH2D* h2_BTaggingEff_csv_loose_Num_udsg = new TH2D("h2_BTaggingEff_csv_loose_Num_udsg", "b-tagging effiency, Numerator, MC, CSV loose, udsg-jets, #eta vs p_{T};p_{T};#eta", Pt_nBins, Pt_bins, Eta_nBins, Eta_bins);

    std::vector<std::string> jetcorr_filenames_25ns_MC_pfL1L2L3;
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("../CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L1FastJet_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("../CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L2Relative_AK4PFchs.txt");
    jetcorr_filenames_25ns_MC_pfL1L2L3.push_back  ("../CORE/Tools/jetcorr/data/run2_25ns/"+jecEraMC+"_MC/"+jecEraMC+"_MC_L3Absolute_AK4PFchs.txt");
    auto jetcorrector = makeJetCorrector(jetcorr_filenames_25ns_MC_pfL1L2L3);

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();

    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);

    tqdm bar;

    while ( (currentFile = (TFile*)fileIter.Next()) ) { 

        if (STOP_REQUESTED) break;

        TString name = ch->GetTitle();
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");
        cms3.Init(tree);

        TString filename(currentFile->GetTitle());

        auto tokens = filename.Tokenize("/");
        auto basename = ((TObjString*)(tokens->At(tokens->GetEntries()-1)))->String().Data();
        if (not quiet) bar.set_label(basename);

        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {
            if (STOP_REQUESTED) break;

            cms3.GetEntry(event);
            nEventsTotal++;

            for (size_t ijet = 0; ijet < pfjets_p4().size(); ++ijet) {
                float cms4pt = pfjets_p4()[ijet].pt();
                float rawpt = cms4pt*tas::pfjets_undoJEC()[ijet];
                float eta = pfjets_p4()[ijet].eta();
                float aeta = fabs(eta);
                jetcorrector->setJetEta(eta); 
                jetcorrector->setJetPt(rawpt); 
                jetcorrector->setJetA(pfjets_area()[ijet]); 
                jetcorrector->setRho(evt_fixgridfastjet_all_rho()); 
                float jec = jetcorrector->getCorrection(); 
                float pt = rawpt*jec;
                if (pt < 20.) continue;
                if (aeta > 2.8) continue;
                int flav = abs(pfjets_hadronFlavour()[ijet]);
                float disc = -1;
                if (is80X) {
                    disc = getbtagvalue("deepFlavourJetTags:probb",ijet) + getbtagvalue("deepFlavourJetTags:probbb",ijet);
                } else {
                    disc = getbtagvalue("pfDeepCSVJetTags:probb",ijet) + getbtagvalue("pfDeepCSVJetTags:probbb",ijet);
                }
                if (flav == 5) { 
                    // std::cout <<  " cms4pt: " << cms4pt <<  " rawpt: " << rawpt <<  " pt: " << pt <<  " jec: " << jec <<  " eta: " << eta <<  std::endl;
                    h2_BTaggingEff_csv_Denom_b->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_TIGHT) h2_BTaggingEff_csv_tight_Num_b->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_MEDIUM) h2_BTaggingEff_csv_med_Num_b->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_LOOSE) h2_BTaggingEff_csv_loose_Num_b->Fill(pt,aeta);
                } else if (flav == 4) {
                    h2_BTaggingEff_csv_Denom_c->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_TIGHT) h2_BTaggingEff_csv_tight_Num_c->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_MEDIUM) h2_BTaggingEff_csv_med_Num_c->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_LOOSE) h2_BTaggingEff_csv_loose_Num_c->Fill(pt,aeta);
                } else if (flav <= 3 || flav == 21) {
                    h2_BTaggingEff_csv_Denom_udsg->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_TIGHT) h2_BTaggingEff_csv_tight_Num_udsg->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_MEDIUM) h2_BTaggingEff_csv_med_Num_udsg->Fill(pt,aeta);
                    if (disc > WP_DEEPCSV_LOOSE) h2_BTaggingEff_csv_loose_Num_udsg->Fill(pt,aeta);
                }
            }

            if (not quiet) bar.progress(nEventsTotal, nEventsChain);

        }//event loop

        delete file;
    }//file loop

    TFile *fout = new TFile(outname, "recreate");
    fout->cd();

    h2_BTaggingEff_csv_Denom_b->Write();
    h2_BTaggingEff_csv_Denom_c->Write();
    h2_BTaggingEff_csv_Denom_udsg->Write();

    h2_BTaggingEff_csv_tight_Num_b->Write();
    h2_BTaggingEff_csv_tight_Num_c->Write();
    h2_BTaggingEff_csv_tight_Num_udsg->Write();

    h2_BTaggingEff_csv_med_Num_b->Write();
    h2_BTaggingEff_csv_med_Num_c->Write();
    h2_BTaggingEff_csv_med_Num_udsg->Write();

    h2_BTaggingEff_csv_loose_Num_b->Write();
    h2_BTaggingEff_csv_loose_Num_c->Write();
    h2_BTaggingEff_csv_loose_Num_udsg->Write();


    fout->Close();

    return 0;

}

