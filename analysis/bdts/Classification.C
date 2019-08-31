#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int Classification( TString myMethodList = "" )
{
   TMVA::Tools::Instance();

   // directory training tree

   TString basedir ="/home/shchauha/2019/Analysis/FCNCAna/analysis/yields/v3.24_training/";
   
   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;      
   // Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start Classification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Here the preparation phase begins

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   TFile * signal = TFile::Open(basedir+"histos_fcnc.root");
   TFile * bg1 = TFile::Open(basedir+"histos_wz.root");
   TFile * bg2 = TFile::Open(basedir+"histos_ttdl.root");
   TFile * bg3 = TFile::Open(basedir+"histos_tth.root");
   TFile * bg4 = TFile::Open(basedir+"histos_ttsl.root");
   TFile * bg5 = TFile::Open(basedir+"histos_ttw.root");	 
   TFile * bg6 = TFile::Open(basedir+"histos_ttz.root");
   
   // Register the training and test trees
   TTree *signalTree     = (TTree*)signal->Get("t");
   TTree *background1     = (TTree*)bg1->Get("t");
   TTree *background2     = (TTree*)bg2->Get("t");
   TTree *background3     = (TTree*)bg3->Get("t");
   TTree *background4     = (TTree*)bg4->Get("t");
   TTree *background5     = (TTree*)bg5->Get("t");
   TTree *background6     = (TTree*)bg6->Get("t");
   //TTree *background7     = (TTree*)bg7->Get("SSTree");

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "BDT_output_trainig_testing.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string

   //actual one
   TMVA::Factory *factory = new TMVA::Factory( "Classification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   //for 1 var, since P doesn't work for 1 var
   //TMVA::Factory *factory = new TMVA::Factory( "Classification", outputFile,
   //                                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]


   dataloader->AddVariable("lep1pt");
   dataloader->AddVariable("lep2pt");
   dataloader->AddVariable("lep1eta");
   dataloader->AddVariable("lep2eta");
   dataloader->AddVariable("njets");
   dataloader->AddVariable("nbtags");
   dataloader->AddVariable("met");
   dataloader->AddVariable("ht");
   dataloader->AddVariable("mll");
   dataloader->AddVariable("nele");
   dataloader->AddVariable("jet1pt");
   dataloader->AddVariable("jet2pt");
   dataloader->AddVariable("btag1pt");
   dataloader->AddVariable("drl1l2");
   dataloader->AddVariable("mindrl1j");
   dataloader->AddVariable("mindrl2j");
   dataloader->AddVariable("mindrl1bt");
   dataloader->AddVariable("mindrl2bt");
   dataloader->AddVariable("dphil1met");
   dataloader->AddVariable("dphil2met");
   dataloader->AddVariable("mt1");
   dataloader->AddVariable("mt2");
   dataloader->AddVariable("dphil1l2");
   dataloader->AddVariable("l1miniiso");
   dataloader->AddVariable("l2miniiso");
   dataloader->AddVariable("l1dxy");
   dataloader->AddVariable("l1dz");
   dataloader->AddVariable("l2dxy");
   dataloader->AddVariable("l2dz");
   dataloader->AddVariable("l1ptratio");
   dataloader->AddVariable("l1ptrel");
   dataloader->AddVariable("l2ptratio");
   dataloader->AddVariable("l2ptrel");
   //dataloader->AddVariable("weight");
   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   // dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   // dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   // dataloader->AddSignalTree    ( signalTree,     signalWeight );
   // dataloader->AddBackgroundTree( background1, backgroundWeight );
   // dataloader->AddBackgroundTree( background2, backgroundWeight );
   // dataloader->AddBackgroundTree( background3, backgroundWeight );
   // dataloader->AddBackgroundTree( background4, backgroundWeight );
   // dataloader->AddBackgroundTree( background5, backgroundWeight );
   // dataloader->AddBackgroundTree( background6, backgroundWeight );

   // To give different trees for training and testing, do as follows:
   //
   dataloader->AddSignalTree( signalTree,     signalWeight, "Training" );
   dataloader->AddSignalTree( signalTree,     signalWeight,  "Test" );
   
   dataloader->AddBackgroundTree( background1, backgroundWeight, "Training" );
   dataloader->AddBackgroundTree( background2, backgroundWeight, "Training" );
   dataloader->AddBackgroundTree( background3, backgroundWeight, "Training" );
   dataloader->AddBackgroundTree( background4, backgroundWeight, "Training" );
   dataloader->AddBackgroundTree( background5, backgroundWeight, "Training" );
   dataloader->AddBackgroundTree( background6, backgroundWeight, "Training" );

   dataloader->AddBackgroundTree( background1, backgroundWeight, "Test" );
   dataloader->AddBackgroundTree( background2, backgroundWeight, "Test" );
   dataloader->AddBackgroundTree( background3, backgroundWeight, "Test" );
   dataloader->AddBackgroundTree( background4, backgroundWeight, "Test" );
   dataloader->AddBackgroundTree( background5, backgroundWeight, "Test" );
   dataloader->AddBackgroundTree( background6, backgroundWeight, "Test" );
            
   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   // ```cpp
   // // --- begin ----------------------------------------------------------
   // std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   // Float_t  treevars[4], weight;
   //
   // // Signal
   // for (UInt_t ivar=0; ivar<4; ivar++) signalTree->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   // for (UInt_t i=0; i<signalTree->GetEntries(); i++) {
   //    signalTree->GetEntry(i);
   //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //    // add training and test events; here: first half is training, second is testing
   //    // note that the weight can also be event-wise
   //    if (i < signalTree->GetEntries()/2.0) dataloader->AddSignalTrainingEvent( vars, signalWeight );
   //    else                              dataloader->AddSignalTestEvent    ( vars, signalWeight );
   // }
   //
   // // Background (has event weights)
   // background->SetBranchAddress( "weight", &weight );
   // for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   // for (UInt_t i=0; i<background->GetEntries(); i++) {
   //    background->GetEntry(i);
   //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //    // add training and test events; here: first half is training, second is testing
   //    // note that the weight can also be event-wise
   //    if (i < background->GetEntries()/2) dataloader->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //    else                                dataloader->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   // }
   // // --- end ------------------------------------------------------------
   // ```
   // End of tree registration

   // Set individual event weights (the variables must exist in the original TTree)
   // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
   // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
   dataloader->SetSignalWeightExpression( "weight" );
   dataloader->SetBackgroundWeightExpression( "weight" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the dataloader how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //
   dataloader->PrepareTrainingAndTestTree( mycuts, "SplitMode=random:!V" );
   //
   // To also specify the number of testing events, use:
   //
   //    dataloader->PrepareTrainingAndTestTree( mycut,
   //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   //dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,"nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );

   // ### Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Other classification techniques deleted from original file were here

   // Boosted Decision Trees

   string config;
   TString name; 

   vector<int> trees = {1000}; //default: 1000
   vector<TString> nodes = {"2.5"}; //default: 2.5%
   vector<int> depth = {2}; //default: 2

   // Gradient Boost
   if (Use["BDTG"]){
        for (int t: trees){
        for (TString n: nodes){
        for (int d: depth){
            config = "!H:!V:NTrees=t:MinNodeSize=n%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=d";
            config.replace(118,1,to_string(d));
            config.replace(27,1,n);
            config.replace(13,1,to_string(t)); 
            name = "BDTG";
            name += to_string(t) + "t" + n + "%n" + to_string(d) + "d";
            factory->BookMethod(dataloader, TMVA::Types::kBDT, name, config);            
        }
        }
        }
    }

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: ClassificationCategory
   //
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   //     factory->OptimizeAllMethods("SigEffAt001","Scan");
   //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> Classification is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return Classification(methodList);
}
