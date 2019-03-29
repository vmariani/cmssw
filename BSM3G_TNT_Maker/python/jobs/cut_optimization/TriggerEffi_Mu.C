/**
This Macro
1. Calls trees, fix cuts and creates ntuples with cut data  

Need to specify
0. See Declare constants
*/
/////
//   To run: root -l SkimNtuple.cc+  
/////
/////
//   Prepare Root and Roofit
/////
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>

using namespace std;

//void filename_()
void TriggerEffi_Mu(){

TChain *a_ = new TChain("BOOM");

a_->Add("/eos/user/v/vmariani/NTuples/DY.root");
//inputFile

int HLT_Mu, HLT_Ele;
std:vector<double>* Muon_pt; Muon_pt=0;
vector<double>* Muon_eta; Muon_eta=0;
vector<double>* Muon_phi; Muon_phi=0;
vector<double>*BoostedJet_pt; BoostedJet_pt=0;
vector<double>*BoostedJet_eta; BoostedJet_eta=0;
vector<double>*BoostedJet_phi; BoostedJet_phi=0;
double numOfHighptEle, numOfVetoEle, numOfHighptMu, numOfLooseMu, numOfBoostedJets;

TBranch *a_HLT_Ele115_CaloIdVT_GsfTrkIdT=a_->GetBranch("HLT_Ele115_CaloIdVT_GsfTrkIdT");
TBranch *a_HLT_Mu50=a_->GetBranch("HLT_Mu50");

TBranch *a_Muon_pt=a_->GetBranch("Muon_pt");

TBranch *a_Muon_eta=a_->GetBranch("Muon_eta");
TBranch *a_Muon_phi=a_->GetBranch("Muon_phi");

TBranch *a_BoostedJet_pt=a_->GetBranch("BoostedJet_pt");
TBranch *a_BoostedJet_eta=a_->GetBranch("BoostedJet_eta");
TBranch *a_BoostedJet_phi=a_->GetBranch("BoostedJet_phi");

TBranch *a_numOfHighptEle=a_->GetBranch("numOfHighptEle");
TBranch *a_numOfHighptMu=a_->GetBranch("numOfHighptMu");
TBranch *a_numOfLooseMu=a_->GetBranch("numOfLooseMu");
TBranch *a_numOfBoostedJets=a_->GetBranch("numOfBoostedJets");
TBranch *a_numOfVetoEle=a_->GetBranch("numOfVetoEle");

a_HLT_Mu50->SetAddress(&HLT_Mu);
a_HLT_Ele115_CaloIdVT_GsfTrkIdT->SetAddress(&HLT_Ele);

a_Muon_pt->SetAddress(&Muon_pt);
a_Muon_eta->SetAddress(&Muon_eta);
a_Muon_phi->SetAddress(&Muon_phi);

a_BoostedJet_pt->SetAddress(&BoostedJet_pt);
a_BoostedJet_eta->SetAddress(&BoostedJet_eta);
a_BoostedJet_phi->SetAddress(&BoostedJet_phi);

a_numOfHighptEle->SetAddress(&numOfHighptEle);
a_numOfHighptMu->SetAddress(&numOfHighptMu);
a_numOfLooseMu->SetAddress(&numOfLooseMu);
a_numOfBoostedJets->SetAddress(&numOfBoostedJets);
a_numOfVetoEle->SetAddress(&numOfVetoEle);

TH1D * ptmu1_050 = new TH1D ("ptmu1_050", "ptmu1_050", 4, 0, 2);
TH1D * ptmu1_5070 = new TH1D ("ptmu1_5070", "ptmu1_5070", 4, 0, 2);
TH1D * ptmu1_70100 = new TH1D ("ptmu1_70100", "ptmu1_70100", 4, 0, 2);
TH1D * ptmu1_100300 = new TH1D ("ptmu1_100300", "ptmu1_100300", 4, 0, 2);
TH1D * ptmu1_300500 = new TH1D ("ptmu1_300500", "ptmu1_300500", 4, 0, 2);
TH1D * ptmu1_500700 = new TH1D ("ptmu1_500700", "ptmu1_500700", 4, 0, 2);
TH1D * ptmu1_7001000 = new TH1D ("ptmu1_7001000", "ptmu1_7001000", 4, 0, 2);
TH1D * ptmu1_10001500 = new TH1D ("ptmu1_10001500", "ptmu1_10001500", 4, 0, 2);
TH1D * ptmu1_15002000 = new TH1D ("ptmu1_15002000", "ptmu1_15002000", 4, 0, 2); 
TH1D * ptmu1_20003000 = new TH1D ("ptmu1_20003000", "ptmu1_20003000", 4, 0, 2);
TH1D * ptmu1_30005000 = new TH1D ("ptmu1_30005000", "ptmu1_30005000", 4, 0, 2);

TH1D * trg_ptmu1_050 = new TH1D ("trg_ptmu1_050", "trg_ptmu1_050", 4, 0, 2);
TH1D * trg_ptmu1_5070 = new TH1D ("trg_ptmu1_5070", "trg_ptmu1_5070", 4, 0, 2);
TH1D * trg_ptmu1_70100 = new TH1D ("trg_ptmu1_70100", "trg_ptmu1_70100", 4, 0, 2);
TH1D * trg_ptmu1_100300 = new TH1D ("trg_ptmu1_100300", "trg_ptmu1_100300", 4, 0, 2);
TH1D * trg_ptmu1_300500 = new TH1D ("trg_ptmu1_300500", "trg_ptmu1_300500", 4, 0, 2);
TH1D * trg_ptmu1_500700 = new TH1D ("trg_ptmu1_500700", "trg_ptmu1_500700", 4, 0, 2);
TH1D * trg_ptmu1_7001000 = new TH1D ("trg_ptmu1_7001000", "trg_ptmu1_7001000", 4, 0, 2);
TH1D * trg_ptmu1_10001500 = new TH1D ("trg_ptmu1_10001500", "trg_ptmu1_10001500", 4, 0, 2);
TH1D * trg_ptmu1_15002000 = new TH1D ("trg_ptmu1_15002000", "trg_ptmu1_15002000", 4, 0, 2);
TH1D * trg_ptmu1_20003000 = new TH1D ("trg_ptmu1_20003000", "trg_ptmu1_20003000", 4, 0, 2);
TH1D * trg_ptmu1_30005000 = new TH1D ("trg_ptmu1_30005000", "trg_ptmu1_30005000", 4, 0, 2);

TH1D * ptmu2_050 = new TH1D ("ptmu2_050", "ptmu2_050", 4, 0, 2);
TH1D * ptmu2_5070 = new TH1D ("ptmu2_5070", "ptmu2_5070", 4, 0, 2);
TH1D * ptmu2_70100 = new TH1D ("ptmu2_70100", "ptmu2_70100", 4, 0, 2);
TH1D * ptmu2_100300 = new TH1D ("ptmu2_100300", "ptmu2_100300", 4, 0, 2);
TH1D * ptmu2_300500 = new TH1D ("ptmu2_300500", "ptmu2_300500", 4, 0, 2);
TH1D * ptmu2_500700 = new TH1D ("ptmu2_500700", "ptmu2_500700", 4, 0, 2);
TH1D * ptmu2_7001000 = new TH1D ("ptmu2_7001000", "ptmu2_7001000", 4, 0, 2);
TH1D * ptmu2_10001500 = new TH1D ("ptmu2_10001500", "ptmu2_10001500", 4, 0, 2);
TH1D * ptmu2_15002000 = new TH1D ("ptmu2_15002000", "ptmu2_15002000", 4, 0, 2);

TH1D * trg_ptmu2_050 = new TH1D ("trg_ptmu2_050", "trg_ptmu2_050", 4, 0, 2);
TH1D * trg_ptmu2_5070 = new TH1D ("trg_ptmu2_5070", "trg_ptmu2_5070", 4, 0, 2);
TH1D * trg_ptmu2_70100 = new TH1D ("trg_ptmu2_70100", "trg_ptmu2_70100", 4, 0, 2);
TH1D * trg_ptmu2_100300 = new TH1D ("trg_ptmu2_100300", "trg_ptmu2_100300", 4, 0, 2);
TH1D * trg_ptmu2_300500 = new TH1D ("trg_ptmu2_300500", "trg_ptmu2_300500", 4, 0, 2);
TH1D * trg_ptmu2_500700 = new TH1D ("trg_ptmu2_500700", "trg_ptmu2_500700", 4, 0, 2);
TH1D * trg_ptmu2_7001000 = new TH1D ("trg_ptmu2_7001000", "trg_ptmu2_7001000", 4, 0, 2);
TH1D * trg_ptmu2_10001500 = new TH1D ("trg_ptmu2_10001500", "trg_ptmu2_10001500", 4, 0, 2);
TH1D * trg_ptmu2_15002000 = new TH1D ("trg_ptmu2_15002000", "trg_ptmu2_15002000", 4, 0, 2);

TH1D *Muon1_pt = new TH1D ("Muon1_pt", "Muon1_pt", 200, 0, 5000);
TH1D *Muon2_pt = new TH1D ("Muon2_pt", "Muon2_pt", 200, 0, 5000);
TH1D *Muon1_phi = new TH1D ("Muon1_phi", "Muon1_phi", 200, -3, 3);
TH1D *Muon2_phi = new TH1D ("Muon2_phi", "Muon2_phi", 200, -3, 3);
TH1D *Muon1_eta = new TH1D ("Muon1_eta", "Muon1_eta", 200, -4, 4);
TH1D *Muon2_eta = new TH1D ("Muon2_eta", "Muon2_eta", 200, -4, 4);

TH1D *Jet_pt = new TH1D ("Jet_pt", "Jet_pt", 200, 0, 5000);
TH1D *Jet_phi = new TH1D ("Jet_phi", "Jet_phi", 200, -3, 3);
TH1D *Jet_eta = new TH1D ("Jet_eta", "Jet_eta", 200, -4, 4);

cout << a_->GetEntries() << endl;
int tot = 0;

for (Int_t i=0;i<a_->GetEntries();i++) {
 a_->GetEntry(i);
 tot = a_->GetEntries();
 if (i%10000 == 0){
 cout << i << " eventi analizzati su " << tot << endl;}
 if (Muon_pt->size() > 0){
  Muon1_pt->Fill(Muon_pt->at(0));
  Muon1_eta->Fill(Muon_eta->at(0));
  Muon1_phi->Fill(Muon_phi->at(0));
  if (Muon_pt->at(0) < 50) ptmu1_050->Fill(1);
  if (Muon_pt->at(0) >= 50 && Muon_pt->at(0) < 70) ptmu1_5070->Fill(1);
  if (Muon_pt->at(0) >= 70 && Muon_pt->at(0) < 100) ptmu1_70100->Fill(1);
  if (Muon_pt->at(0) >= 100 && Muon_pt->at(0) < 300) ptmu1_100300->Fill(1);
  if (Muon_pt->at(0) >= 300 && Muon_pt->at(0) < 500) ptmu1_300500->Fill(1);
  if (Muon_pt->at(0) >= 500 && Muon_pt->at(0) < 700) ptmu1_500700->Fill(1);
  if (Muon_pt->at(0) >= 700 && Muon_pt->at(0) < 1000) ptmu1_7001000->Fill(1);
  if (Muon_pt->at(0) >= 1000 && Muon_pt->at(0) < 1500) ptmu1_10001500->Fill(1);
  if (Muon_pt->at(0) >= 1500 && Muon_pt->at(0) < 2000) ptmu1_15002000->Fill(1);
  if (Muon_pt->at(0) >= 2000 && Muon_pt->at(0) < 3000) ptmu1_20003000->Fill(1);
  if (Muon_pt->at(0) >= 3000 && Muon_pt->at(0) < 5000) ptmu1_30005000->Fill(1);  
  if (Muon_pt->size() > 1){
   Muon2_pt->Fill(Muon_pt->at(1));
   Muon2_eta->Fill(Muon_eta->at(1));
   Muon2_phi->Fill(Muon_phi->at(1));
   if (BoostedJet_pt->size() > 0){
    Jet_pt->Fill(BoostedJet_pt->at(0));
    Jet_eta->Fill(BoostedJet_eta->at(0));
    Jet_phi->Fill(BoostedJet_phi->at(0));
   }
   if (Muon_pt->at(1) < 50) ptmu2_050->Fill(1);
   if (Muon_pt->at(1) >= 50 && Muon_pt->at(1) < 70) ptmu2_5070->Fill(1);
   if (Muon_pt->at(1) >= 70 && Muon_pt->at(1) < 100) ptmu2_70100->Fill(1);
   if (Muon_pt->at(1) >= 100 && Muon_pt->at(1) < 300) ptmu2_100300->Fill(1);
   if (Muon_pt->at(1) >= 300 && Muon_pt->at(1) < 500) ptmu2_300500->Fill(1);
   if (Muon_pt->at(1) >= 500 && Muon_pt->at(1) < 700) ptmu2_500700->Fill(1);
   if (Muon_pt->at(1) >= 700 && Muon_pt->at(1) < 1000) ptmu2_7001000->Fill(1);
   if (Muon_pt->at(1) >= 1000 && Muon_pt->at(1) < 1500) ptmu2_10001500->Fill(1);
   if (Muon_pt->at(1) >= 1500 && Muon_pt->at(1) < 2000) ptmu2_15002000->Fill(1); 
   }
 
  if (HLT_Mu > 0){
   if (Muon_pt->at(0) < 50)trg_ptmu1_050->Fill(1);
   if (Muon_pt->at(0) >= 50 && Muon_pt->at(0) < 70)trg_ptmu1_5070->Fill(1);
   if (Muon_pt->at(0) >= 70 && Muon_pt->at(0) < 100)trg_ptmu1_70100->Fill(1);
   if (Muon_pt->at(0) >= 100 && Muon_pt->at(0) < 300)trg_ptmu1_100300->Fill(1);
   if (Muon_pt->at(0) >= 300 && Muon_pt->at(0) < 500)trg_ptmu1_300500->Fill(1);
   if (Muon_pt->at(0) >= 500 && Muon_pt->at(0) < 700)trg_ptmu1_500700->Fill(1);
   if (Muon_pt->at(0) >= 700 && Muon_pt->at(0) < 1000)trg_ptmu1_7001000->Fill(1);
   if (Muon_pt->at(0) >= 1000 && Muon_pt->at(0) < 1500)trg_ptmu1_10001500->Fill(1);   
   if (Muon_pt->at(0) >= 1500 && Muon_pt->at(0) < 2000)trg_ptmu1_15002000->Fill(1);
   if (Muon_pt->at(0) >= 2000 && Muon_pt->at(0) < 3000)trg_ptmu1_20003000->Fill(1);
   if (Muon_pt->at(0) >= 3000 && Muon_pt->at(0) < 5000)trg_ptmu1_30005000->Fill(1);
  
   if (Muon_pt->size() > 1){
    if (Muon_pt->at(1) < 50)trg_ptmu2_050->Fill(1);
    if (Muon_pt->at(1) >= 50 && Muon_pt->at(1) < 70)trg_ptmu2_5070->Fill(1);
    if (Muon_pt->at(1) >= 70 && Muon_pt->at(1) < 100)trg_ptmu2_70100->Fill(1);
    if (Muon_pt->at(1) >= 100 && Muon_pt->at(1) < 300)trg_ptmu2_100300->Fill(1);
    if (Muon_pt->at(1) >= 300 && Muon_pt->at(1) < 500)trg_ptmu2_300500->Fill(1);
    if (Muon_pt->at(1) >= 500 && Muon_pt->at(1) < 700)trg_ptmu2_500700->Fill(1);
    if (Muon_pt->at(1) >= 700 && Muon_pt->at(1) < 1000)trg_ptmu2_7001000->Fill(1);
    if (Muon_pt->at(1) >= 1000 && Muon_pt->at(1) < 1500)trg_ptmu2_10001500->Fill(1); 
    if (Muon_pt->at(1) >= 1500 && Muon_pt->at(1) < 2000)trg_ptmu2_15002000->Fill(1);
   }
  }
 }
}

TFile *f = new TFile("Trigger_effi_mu_DY.root", "RECREATE");
Muon1_pt->Write();
Muon2_pt->Write();
Muon1_eta->Write();
Muon2_eta->Write();
Muon1_phi->Write();
Muon2_phi->Write();
Jet_pt->Write();
Jet_eta->Write();
Jet_phi->Write();
ptmu1_050->Write();
ptmu1_5070->Write();
ptmu1_70100->Write();
ptmu1_100300->Write();
ptmu1_300500->Write();
ptmu1_500700->Write();
ptmu1_7001000->Write();
ptmu1_10001500->Write();
ptmu1_15002000->Write();
ptmu1_20003000->Write();
ptmu1_30005000->Write();

trg_ptmu1_050->Write();
trg_ptmu1_5070->Write();
trg_ptmu1_70100->Write();
trg_ptmu1_100300->Write();
trg_ptmu1_300500->Write();
trg_ptmu1_500700->Write();
trg_ptmu1_7001000->Write();
trg_ptmu1_10001500->Write();
trg_ptmu1_15002000->Write();
trg_ptmu1_20003000->Write();
trg_ptmu1_30005000->Write();

ptmu2_050->Write();
ptmu2_5070->Write();
ptmu2_70100->Write();
ptmu2_100300->Write();
ptmu2_300500->Write();
ptmu2_500700->Write();
ptmu2_7001000->Write();
ptmu2_10001500->Write();
ptmu2_15002000->Write();

trg_ptmu2_050->Write();
trg_ptmu2_5070->Write();
trg_ptmu2_70100->Write();
trg_ptmu2_100300->Write();
trg_ptmu2_300500->Write();
trg_ptmu2_500700->Write();
trg_ptmu2_7001000->Write();
trg_ptmu2_10001500->Write();
trg_ptmu2_15002000->Write();


f->Write();
f->Close();

}
 
