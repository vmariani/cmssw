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
void Analisi_Mu_BDT(){

TChain *a_ = new TChain("BOOM");

a_->Add("/eos/user/v/vmariani/NTuples/HN/ST_tot.root");

int HLT_Mu;
std:vector<double>* Muon_pt; Muon_pt=0;
vector<double>* Muon_eta; Muon_eta=0;
vector<double>* Muon_phi; Muon_phi=0;
vector<double>*BoostedJet_pt; BoostedJet_pt=0;
vector<double>*BoostedJet_eta; BoostedJet_eta=0;
vector<double>*BoostedJet_phi; BoostedJet_phi=0;
double M_leplep;
double numOfHighptEle, numOfVetoEle, numOfHighptMu, numOfLooseMu, numOfBoostedJets;
double lepsf_evt, lumi_wgt, musf_trigger_mu1 ,musf_ID_mu1, musf_iso_mu1, musf_tot_mu1, musf_trigger_mu2, musf_ID_mu2, musf_iso_mu2, musf_tot_mu2;

double Muon1_pt, Muon1_eta, Muon1_phi, Muon2_pt, Muon2_eta, Muon2_phi, Jet_pt, Jet_eta ,Jet_phi, M_ll, weight;

TTree *Tree = new TTree("Analysis","");

Tree->Branch("Muon1_pt",&Muon1_pt, "Muon1_pt/D");
Tree->Branch("Muon1_eta",&Muon1_eta, "Muon1_eta/D");
Tree->Branch("Muon1_phi",&Muon1_phi, "Muon1_phi/D");

Tree->Branch("Muon2_pt",&Muon2_pt, "Muon2_pt/D");
Tree->Branch("Muon2_eta",&Muon2_eta, "Muon2_eta/D");
Tree->Branch("Muon2_phi",&Muon2_phi, "Muon2_phi/D");

Tree->Branch("Jet_pt",&Jet_pt, "Jet_pt/D");
Tree->Branch("Jet_eta",&Jet_eta, "Jet_eta/D");
Tree->Branch("Jet_phi",&Jet_phi, "Jet_phi/D");

Tree->Branch("M_ll", &M_ll, "M_ll/D");
Tree->Branch("weight", &weight, "weight/D");


TBranch *a_HLT_Mu50=a_->GetBranch("HLT_Mu50");

TBranch *a_Muon_pt=a_->GetBranch("Muon_pt");
TBranch *a_Muon_eta=a_->GetBranch("Muon_eta");
TBranch *a_Muon_phi=a_->GetBranch("Muon_phi");

TBranch *a_M_leplep=a_->GetBranch("M_leplep");

TBranch *a_lepsf_evt=a_->GetBranch("lepsf_evt");
TBranch *a_musf_trigger_mu1=a_->GetBranch("musf_trigger_mu1");
TBranch *a_musf_ID_mu1=a_->GetBranch("musf_ID_mu1");
TBranch *a_musf_iso_mu1=a_->GetBranch("musf_iso_mu1");
TBranch *a_musf_tot_mu1=a_->GetBranch("musf_tot_mu1");
TBranch *a_musf_trigger_mu2=a_->GetBranch("musf_trigger_mu2");
TBranch *a_musf_ID_mu2=a_->GetBranch("musf_ID_mu2");
TBranch *a_musf_iso_mu2=a_->GetBranch("musf_iso_mu2");
TBranch *a_musf_tot_mu2=a_->GetBranch("musf_tot_mu2");
TBranch *a_lumi_wgt=a_->GetBranch("lumi_wgt");

TBranch *a_BoostedJet_pt=a_->GetBranch("BoostedJet_pt");
TBranch *a_BoostedJet_eta=a_->GetBranch("BoostedJet_eta");
TBranch *a_BoostedJet_phi=a_->GetBranch("BoostedJet_phi");

TBranch *a_numOfHighptEle=a_->GetBranch("numOfHighptEle");
TBranch *a_numOfHighptMu=a_->GetBranch("numOfHighptMu");
TBranch *a_numOfLooseMu=a_->GetBranch("numOfLooseMu");
TBranch *a_numOfBoostedJets=a_->GetBranch("numOfBoostedJets");
TBranch *a_numOfVetoEle=a_->GetBranch("numOfVetoEle");

a_HLT_Mu50->SetAddress(&HLT_Mu);

a_Muon_pt->SetAddress(&Muon_pt);
a_Muon_eta->SetAddress(&Muon_eta);
a_Muon_phi->SetAddress(&Muon_phi);

a_M_leplep->SetAddress(&M_leplep);

a_lepsf_evt->SetAddress(&lepsf_evt);
a_musf_trigger_mu1->SetAddress(&musf_trigger_mu1);
a_musf_ID_mu1->SetAddress(&musf_ID_mu1);
a_musf_iso_mu1->SetAddress(&musf_iso_mu1);
a_musf_tot_mu1->SetAddress(&musf_tot_mu1);
a_musf_trigger_mu2->SetAddress(&musf_trigger_mu2);
a_musf_ID_mu2->SetAddress(&musf_ID_mu2);
a_musf_iso_mu2->SetAddress(&musf_iso_mu2);
a_musf_tot_mu2->SetAddress(&musf_tot_mu2);
a_lumi_wgt->SetAddress(&lumi_wgt);

a_BoostedJet_pt->SetAddress(&BoostedJet_pt);
a_BoostedJet_eta->SetAddress(&BoostedJet_eta);
a_BoostedJet_phi->SetAddress(&BoostedJet_phi);

a_numOfHighptEle->SetAddress(&numOfHighptEle);
a_numOfHighptMu->SetAddress(&numOfHighptMu);
a_numOfLooseMu->SetAddress(&numOfLooseMu);
a_numOfBoostedJets->SetAddress(&numOfBoostedJets);
a_numOfVetoEle->SetAddress(&numOfVetoEle);

cout << a_->GetEntries() << endl;
int tot = 0;

double wg = 0;
double lumi = 41529;
for (Int_t i=0;i<a_->GetEntries();i++){ 
 a_->GetEntry(i);
 tot = a_->GetEntries();
 if (i%100000 == 0)cout << i << " eventi analizzati su " << tot << endl;
 if (Muon_pt->size() > 1 && BoostedJet_pt->size() > 0){

  if (numOfHighptMu==2 && numOfVetoEle==0 && numOfBoostedJets>= 1 && HLT_Mu==1 && Muon_pt->at(0) > 70 && Muon_pt->at(1) > 50 && fabs(Muon_eta->at(0)) < 2.4 && fabs(Muon_eta->at(1)) < 2.4 && BoostedJet_pt->at(0) > 190 && M_leplep > 200){
   wg = lumi*lumi_wgt*lepsf_evt;
   //printf("wg:%f lepsf_evt:%f lumi_wgt:%f lumi:%f \n", wg, lepsf_evt, lumi_wgt, lumi);

   weight = wg;
   Muon1_pt = Muon_pt->at(0);
   Muon1_eta = Muon_eta->at(0);
   Muon1_phi = Muon_phi->at(0);

   Muon2_pt = Muon_pt->at(1);
   Muon2_eta = Muon_eta->at(1);
   Muon2_phi = Muon_phi->at(1);
 
   Jet_pt= BoostedJet_pt->at(0);
   Jet_eta= BoostedJet_eta->at(0);
   Jet_phi= BoostedJet_phi->at(0);

   M_ll = M_leplep;

   Tree->Fill(); 
  }
 }
}

TFile *f = new TFile("Muon_BDT_ST_tot.root", "RECREATE");
Tree->Write();

f->Write();
f->Close();

}
 
