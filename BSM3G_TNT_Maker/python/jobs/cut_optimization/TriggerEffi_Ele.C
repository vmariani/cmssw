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
void TriggerEffi_Ele(){

TChain *a_ = new TChain("BOOM");

a_->Add("/eos/user/v/vmariani/NTuples/Sign_elel_L7_M7000.root");
//inputFile

int HLT_Mu, HLT_Ele;
std:vector<double>* patElectron_pt; patElectron_pt=0;
vector<double>* patElectron_eta; patElectron_eta=0;
vector<double>* patElectron_phi; patElectron_phi=0;
vector<double>*BoostedJet_pt; BoostedJet_pt=0;
vector<double>*BoostedJet_eta; BoostedJet_eta=0;
vector<double>*BoostedJet_phi; BoostedJet_phi=0;
double numOfHighptEle, numOfVetoEle, numOfHighptMu, numOfLooseMu, numOfBoostedJets;

TBranch *a_HLT_Ele115_CaloIdVT_GsfTrkIdT=a_->GetBranch("HLT_Ele115_CaloIdVT_GsfTrkIdT");
TBranch *a_HLT_Mu50=a_->GetBranch("HLT_Mu50");

TBranch *a_patElectron_pt=a_->GetBranch("patElectron_pt");
TBranch *a_patElectron_eta=a_->GetBranch("patElectron_eta");
TBranch *a_patElectron_phi=a_->GetBranch("patElectron_phi");

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

a_patElectron_pt->SetAddress(&patElectron_pt);
a_patElectron_eta->SetAddress(&patElectron_eta);
a_patElectron_phi->SetAddress(&patElectron_phi);

a_BoostedJet_pt->SetAddress(&BoostedJet_pt);
a_BoostedJet_eta->SetAddress(&BoostedJet_eta);
a_BoostedJet_phi->SetAddress(&BoostedJet_phi);

a_numOfHighptEle->SetAddress(&numOfHighptEle);
a_numOfHighptMu->SetAddress(&numOfHighptMu);
a_numOfLooseMu->SetAddress(&numOfLooseMu);
a_numOfBoostedJets->SetAddress(&numOfBoostedJets);
a_numOfVetoEle->SetAddress(&numOfVetoEle);

TH1D * ptele1_050 = new TH1D ("ptele1_050", "ptele1_050", 4, 0, 2);
TH1D * ptele1_5070 = new TH1D ("ptele1_5070", "ptele1_5070", 4, 0, 2);
TH1D * ptele1_70100 = new TH1D ("ptele1_70100", "ptele1_70100", 4, 0, 2);
TH1D * ptele1_100300 = new TH1D ("ptele1_100300", "ptele1_100300", 4, 0, 2);
TH1D * ptele1_150500 = new TH1D ("ptele1_150500", "ptele1_150500", 4, 0, 2);
TH1D * ptele1_200300 = new TH1D ("ptele1_200300", "ptele1_200300", 4, 0, 2);
TH1D * ptele1_300500 = new TH1D ("ptele1_300500", "ptele1_300500", 4, 0, 2);
TH1D * ptele1_400500 = new TH1D ("ptele1_400500", "ptele1_400500", 4, 0, 2);
TH1D * ptele1_500700 = new TH1D ("ptele1_500700", "ptele1_500700", 4, 0, 2);
TH1D * ptele1_7001000 = new TH1D ("ptele1_7001000", "ptele1_7001000", 4, 0, 2);
TH1D * ptele1_10001500 = new TH1D ("ptele1_10001500", "ptele1_10001500", 4, 0, 2);
TH1D * ptele1_15005000 = new TH1D ("ptele1_15005000", "ptele1_15005000", 4, 0, 2); 
TH1D * ptele1_20003000 = new TH1D ("ptele1_20003000", "ptele1_20003000", 4, 0, 2);
TH1D * ptele1_30005000 = new TH1D ("ptele1_30005000", "ptele1_30005000", 4, 0, 2);

TH1D * trg_ptele1_050 = new TH1D ("trg_ptele1_050", "trg_ptele1_050", 4, 0, 2);
TH1D * trg_ptele1_5070 = new TH1D ("trg_ptele1_5070", "trg_ptele1_5070", 4, 0, 2);
TH1D * trg_ptele1_70100 = new TH1D ("trg_ptele1_70100", "trg_ptele1_70100", 4, 0, 2);
TH1D * trg_ptele1_100300 = new TH1D ("trg_ptele1_100300", "trg_ptele1_100300", 4, 0, 2);
TH1D * trg_ptele1_150500 = new TH1D ("trg_ptele1_150500", "trg_ptele1_150500", 4, 0, 2);
TH1D * trg_ptele1_200300 = new TH1D ("trg_ptele1_200300", "trg_ptele1_200300", 4, 0, 2);
TH1D * trg_ptele1_300500 = new TH1D ("trg_ptele1_300500", "trg_ptele1_300500", 4, 0, 2);
TH1D * trg_ptele1_400500 = new TH1D ("trg_ptele1_400500", "trg_ptele1_400500", 4, 0, 2);
TH1D * trg_ptele1_500700 = new TH1D ("trg_ptele1_500700", "trg_ptele1_500700", 4, 0, 2);
TH1D * trg_ptele1_7001000 = new TH1D ("trg_ptele1_7001000", "trg_ptele1_7001000", 4, 0, 2);
TH1D * trg_ptele1_10001500 = new TH1D ("trg_ptele1_10001500", "trg_ptele1_10001500", 4, 0, 2);
TH1D * trg_ptele1_15005000 = new TH1D ("trg_ptele1_15005000", "trg_ptele1_15005000", 4, 0, 2);
TH1D * trg_ptele1_20003000 = new TH1D ("trg_ptele1_20003000", "trg_ptele1_20003000", 4, 0, 2);
TH1D * trg_ptele1_30005000 = new TH1D ("trg_ptele1_30005000", "trg_ptele1_30005000", 4, 0, 2);

TH1D * ptele2_050 = new TH1D ("ptele2_050", "ptele2_050", 4, 0, 2);
TH1D * ptele2_5070 = new TH1D ("ptele2_5070", "ptele2_5070", 4, 0, 2);
TH1D * ptele2_70100 = new TH1D ("ptele2_70100", "ptele2_70100", 4, 0, 2);
TH1D * ptele2_100300 = new TH1D ("ptele2_100300", "ptele2_100300", 4, 0, 2);
TH1D * ptele2_150500 = new TH1D ("ptele2_150500", "ptele2_150500", 4, 0, 2);
TH1D * ptele2_200300 = new TH1D ("ptele2_200300", "ptele2_200300", 4, 0, 2);
TH1D * ptele2_300500 = new TH1D ("ptele2_300500", "ptele2_300500", 4, 0, 2);
TH1D * ptele2_400500 = new TH1D ("ptele2_400500", "ptele2_400500", 4, 0, 2);
TH1D * ptele2_500700 = new TH1D ("ptele2_500700", "ptele2_500700", 4, 0, 2);
TH1D * ptele2_7001000 = new TH1D ("ptele2_7001000", "ptele2_7001000", 4, 0, 2);
TH1D * ptele2_10001500 = new TH1D ("ptele2_10001500", "ptele2_10001500", 4, 0, 2);
TH1D * ptele2_15005000 = new TH1D ("ptele2_15005000", "ptele2_15005000", 4, 0, 2);

TH1D * trg_ptele2_050 = new TH1D ("trg_ptele2_050", "trg_ptele2_050", 4, 0, 2);
TH1D * trg_ptele2_5070 = new TH1D ("trg_ptele2_5070", "trg_ptele2_5070", 4, 0, 2);
TH1D * trg_ptele2_70100 = new TH1D ("trg_ptele2_70100", "trg_ptele2_70100", 4, 0, 2);
TH1D * trg_ptele2_100300 = new TH1D ("trg_ptele2_100300", "trg_ptele2_100300", 4, 0, 2);
TH1D * trg_ptele2_150500 = new TH1D ("trg_ptele2_150500", "trg_ptele2_150500", 4, 0, 2);
TH1D * trg_ptele2_200300 = new TH1D ("trg_ptele2_200300", "trg_ptele2_200300", 4, 0, 2);
TH1D * trg_ptele2_300500 = new TH1D ("trg_ptele2_300500", "trg_ptele2_300500", 4, 0, 2);
TH1D * trg_ptele2_400500 = new TH1D ("trg_ptele2_400500", "trg_ptele2_400500", 4, 0, 2);
TH1D * trg_ptele2_500700 = new TH1D ("trg_ptele2_500700", "trg_ptele2_500700", 4, 0, 2);
TH1D * trg_ptele2_7001000 = new TH1D ("trg_ptele2_7001000", "trg_ptele2_7001000", 4, 0, 2);
TH1D * trg_ptele2_10001500 = new TH1D ("trg_ptele2_10001500", "trg_ptele2_10001500", 4, 0, 2);
TH1D * trg_ptele2_15005000 = new TH1D ("trg_ptele2_15005000", "trg_ptele2_15005000", 4, 0, 2);

TH1D *Ele1_pt = new TH1D ("Ele1_pt", "Ele1_pt", 200, 0, 5000);
TH1D *Ele2_pt = new TH1D ("Ele2_pt", "Ele2_pt", 200, 0, 5000);
TH1D *Ele1_phi = new TH1D ("Ele1_phi", "Ele1_phi", 200, -3, 3);
TH1D *Ele2_phi = new TH1D ("Ele2_phi", "Ele2_phi", 200, -3, 3);
TH1D *Ele1_eta = new TH1D ("Ele1_eta", "Ele1_eta", 200, -4, 4);
TH1D *Ele2_eta = new TH1D ("Ele2_eta", "Ele2_eta", 200, -4, 4);

TH1D *Jet_pt = new TH1D ("Jet_pt", "Jet_pt", 200, 0, 5000);
TH1D *Jet_phi = new TH1D ("Jet_phi", "Jet_phi", 200, -3, 3);
TH1D *Jet_eta = new TH1D ("Jet_eta", "Jet_eta", 200, -4, 4);

cout << a_->GetEntries() << endl;
int tot=0;

for (Int_t i=0;i<a_->GetEntries();i++) {
 a_->GetEntry(i);
 tot = a_->GetEntries();
 if (i%10000 == 0){
 cout << i << " eventi analizzati su " << tot << endl;}
 if (patElectron_pt->size() > 0){
  if (patElectron_pt->at(0) < 50) ptele1_050->Fill(1);
  if (patElectron_pt->at(0) >= 50 && patElectron_pt->at(0) < 70) ptele1_5070->Fill(1); 
  if (patElectron_pt->at(0) >= 70 && patElectron_pt->at(0) < 100) ptele1_70100->Fill(1);
  if (patElectron_pt->at(0) >= 100 && patElectron_pt->at(0) < 300) ptele1_100300->Fill(1);
  if (patElectron_pt->at(0) >= 300 && patElectron_pt->at(0) < 500) ptele1_300500->Fill(1);
  if (patElectron_pt->at(0) >= 500 && patElectron_pt->at(0) < 700) ptele1_500700->Fill(1);
  if (patElectron_pt->at(0) >= 700 && patElectron_pt->at(0) < 1000) ptele1_7001000->Fill(1);
  if (patElectron_pt->at(0) >= 1000 && patElectron_pt->at(0) < 1500) ptele1_10001500->Fill(1);
  if (patElectron_pt->at(0) >= 1500 && patElectron_pt->at(0) < 2000) ptele1_15005000->Fill(1);
  if (patElectron_pt->at(0) >= 2000 && patElectron_pt->at(0) < 3000) ptele1_20003000->Fill(1);
  if (patElectron_pt->at(0) >= 3000 && patElectron_pt->at(0) < 5000) ptele1_30005000->Fill(1);

  if (patElectron_pt->size() > 1){
   Ele1_pt->Fill(patElectron_pt->at(0));
   Ele1_eta->Fill(patElectron_eta->at(0));
   Ele1_phi->Fill(patElectron_phi->at(0));
   if (patElectron_pt->at(1) < 50) ptele2_050->Fill(1);
   if (patElectron_pt->at(1) >= 50 && patElectron_pt->at(1) < 70) ptele2_5070->Fill(1);
   if (patElectron_pt->at(1) >= 70 && patElectron_pt->at(1) < 100) ptele2_70100->Fill(1);
   if (patElectron_pt->at(1) >= 100 && patElectron_pt->at(1) < 300) ptele2_100300->Fill(1);
   if (patElectron_pt->at(1) >= 300 && patElectron_pt->at(1) < 500) ptele2_300500->Fill(1);
   if (patElectron_pt->at(1) >= 500 && patElectron_pt->at(1) < 700) ptele2_500700->Fill(1);
   if (patElectron_pt->at(1) >= 700 && patElectron_pt->at(1) < 1000) ptele2_7001000->Fill(1);
   if (patElectron_pt->at(1) >= 1000 && patElectron_pt->at(1) < 1500) ptele2_10001500->Fill(1);
   if (patElectron_pt->at(1) >= 1500 && patElectron_pt->at(1) < 2000) ptele2_15005000->Fill(1); 
  }  
 
  if (HLT_Ele > 0){
   if (patElectron_pt->at(0) < 50)trg_ptele1_050->Fill(1);
   if (patElectron_pt->at(0) >= 50 && patElectron_pt->at(0) < 70)trg_ptele1_5070->Fill(1);
   if (patElectron_pt->at(0) >= 70 && patElectron_pt->at(0) < 100)trg_ptele1_70100->Fill(1);
   if (patElectron_pt->at(0) >= 100 && patElectron_pt->at(0) < 300)trg_ptele1_100300->Fill(1);
   if (patElectron_pt->at(0) >= 300 && patElectron_pt->at(0) < 500)trg_ptele1_300500->Fill(1);
   if (patElectron_pt->at(0) >= 500 && patElectron_pt->at(0) < 700)trg_ptele1_500700->Fill(1);
   if (patElectron_pt->at(0) >= 700 && patElectron_pt->at(0) < 1000)trg_ptele1_7001000->Fill(1);
   if (patElectron_pt->at(0) >= 1000 && patElectron_pt->at(0) < 1500)trg_ptele1_10001500->Fill(1);   
   if (patElectron_pt->at(0) >= 1500 && patElectron_pt->at(0) < 2000)trg_ptele1_15005000->Fill(1);
   if (patElectron_pt->at(0) >= 2000 && patElectron_pt->at(0) < 3000)trg_ptele1_20003000->Fill(1);
   if (patElectron_pt->at(0) >= 3000 && patElectron_pt->at(0) < 5000)trg_ptele1_30005000->Fill(1);
 
   if (patElectron_pt->size() > 1){ 
    Ele2_pt->Fill(patElectron_pt->at(1));
    Ele2_eta->Fill(patElectron_eta->at(1));
    Ele2_phi->Fill(patElectron_phi->at(1));
    if (BoostedJet_pt->size() > 0){
     Jet_pt->Fill(BoostedJet_pt->at(0));
     Jet_eta->Fill(BoostedJet_eta->at(0));
     Jet_phi->Fill(BoostedJet_phi->at(0));
    }
    if (patElectron_pt->at(1) < 50)trg_ptele2_050->Fill(1);
    if (patElectron_pt->at(1) >= 50 && patElectron_pt->at(1) < 70)trg_ptele2_5070->Fill(1);
    if (patElectron_pt->at(1) >= 70 && patElectron_pt->at(1) < 100)trg_ptele2_70100->Fill(1);
    if (patElectron_pt->at(1) >= 100 && patElectron_pt->at(1) < 300)trg_ptele2_100300->Fill(1);
    if (patElectron_pt->at(1) >= 300 && patElectron_pt->at(1) < 500)trg_ptele2_300500->Fill(1);
    if (patElectron_pt->at(1) >= 500 && patElectron_pt->at(1) < 700)trg_ptele2_500700->Fill(1);
    if (patElectron_pt->at(1) >= 700 && patElectron_pt->at(1) < 1000)trg_ptele2_7001000->Fill(1);
    if (patElectron_pt->at(1) >= 1000 && patElectron_pt->at(1) < 1500)trg_ptele2_10001500->Fill(1); 
    if (patElectron_pt->at(1) >= 1500 && patElectron_pt->at(1) < 2000)trg_ptele2_15005000->Fill(1);
   }
  }
 }
}

TFile *f = new TFile("Trigger_effi_ele_Sign_elel_L7_M7000.root", "RECREATE");
Ele1_pt->Write();
Ele2_pt->Write();
Ele1_eta->Write();
Ele2_eta->Write();
Ele1_phi->Write();
Ele2_phi->Write();
Jet_pt->Write();
Jet_eta->Write();
Jet_phi->Write();
ptele1_050->Write();
ptele1_5070->Write();
ptele1_70100->Write();
ptele1_100300->Write();
ptele1_300500->Write();
ptele1_500700->Write();
ptele1_7001000->Write();
ptele1_10001500->Write();
ptele1_15005000->Write();
ptele1_20003000->Write();
ptele1_30005000->Write();

trg_ptele1_050->Write();
trg_ptele1_5070->Write();
trg_ptele1_70100->Write();
trg_ptele1_100300->Write();
trg_ptele1_300500->Write();
trg_ptele1_500700->Write();
trg_ptele1_7001000->Write();
trg_ptele1_10001500->Write();
trg_ptele1_15005000->Write();
trg_ptele1_20003000->Write();
trg_ptele1_30005000->Write();

ptele2_050->Write();
ptele2_5070->Write();
ptele2_70100->Write();
ptele2_100300->Write();
ptele2_300500->Write();
ptele2_500700->Write();
ptele2_7001000->Write();
ptele2_10001500->Write();
ptele2_15005000->Write();

trg_ptele2_050->Write();
trg_ptele2_5070->Write();
trg_ptele2_70100->Write();
trg_ptele2_100300->Write();
trg_ptele2_300500->Write();
trg_ptele2_500700->Write();
trg_ptele2_7001000->Write();
trg_ptele2_10001500->Write();
trg_ptele2_15005000->Write();


f->Write();
f->Close();

}
 
