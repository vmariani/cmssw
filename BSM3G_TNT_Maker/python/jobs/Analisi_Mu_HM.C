# define M_PI           3.14159265358979323846  /* pi */
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
void Analisi_Mu_HM(){

TChain *a_ = new TChain("BOOM");

a_->Add("/eos/user/v/vmariani/NTuples/HN/Sign_mumu_L11_M7000.root");

std:vector<double>* Muon_pt; Muon_pt=0;
vector<double>* Muon_eta; Muon_eta=0;
vector<double>* Muon_phi; Muon_phi=0;
vector<double>*BoostedJet_pt; BoostedJet_pt=0;
vector<double>*BoostedJet_eta; BoostedJet_eta=0;
vector<double>*BoostedJet_phi; BoostedJet_phi=0;
vector<double>*Jet_pt;Jet_pt=0;
vector<double>*Jet_eta;Jet_eta=0;
vector<double>*Jet_phi;Jet_phi=0;
int HLT_Mu;
double M_leplep;
double M_leplepBjet;
double numOfHighptEle, numOfVetoEle, numOfHighptMu, numOfLooseMu, numOfBoostedJets;
double lepsf_evt, lumi_wgt, musf_trigger_mu1 ,musf_ID_mu1, musf_iso_mu1, musf_tot_mu1, musf_trigger_mu2, musf_ID_mu2, musf_iso_mu2, musf_tot_mu2;

TBranch *a_HLT_Mu50=a_->GetBranch("HLT_Mu50");
TBranch *a_HLT_OldMu100=a_->GetBranch("HLT_OldMu100");

TBranch *a_Muon_pt=a_->GetBranch("Muon_pt");
TBranch *a_Muon_eta=a_->GetBranch("Muon_eta");
TBranch *a_Muon_phi=a_->GetBranch("Muon_phi");

TBranch *a_M_leplep=a_->GetBranch("M_leplep");
TBranch *a_M_leplepBjet=a_->GetBranch("M_leplepBjet");

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

TBranch *a_Jet_pt=a_->GetBranch("Jet_pt");
TBranch *a_Jet_eta=a_->GetBranch("Jet_eta");
TBranch *a_Jet_phi=a_->GetBranch("Jet_phi");

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
a_M_leplepBjet->SetAddress(&M_leplepBjet);

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
a_Jet_pt->SetAddress(&Jet_pt);
a_Jet_eta->SetAddress(&Jet_eta);
a_Jet_phi->SetAddress(&Jet_phi);

a_numOfHighptEle->SetAddress(&numOfHighptEle);
a_numOfHighptMu->SetAddress(&numOfHighptMu);
a_numOfLooseMu->SetAddress(&numOfLooseMu);
a_numOfBoostedJets->SetAddress(&numOfBoostedJets);
a_numOfVetoEle->SetAddress(&numOfVetoEle);

Float_t bins2[] = {0,200,400,600,800,1000,1400,2000,3500,10000};
Int_t  binnum2 = 9;

//TH1F* h = new TH1F("n","t", binnum, bins);
//
TH1D * M_lljet = new TH1D ("M_lljet", "M_lljet", binnum2, bins2);
//TH1D * DSetagen = new TH1D ("DSetagen", "DSetagen", binnum2, bins2);

TH1D *Jet_separation = new TH1D ("Jet_separation", "Jet_separation", 100, 0, 5);
TH2D *Jet_separation_vs_pt = new TH2D ("Jet_separation_vs_pt", "Jet_separation_vs_pt", 50, 0, 5, 100, 0, 10000);

TH1D *num = new TH1D ("num", "num", 4, 0, 2);
//TH1D *M_lljet = new TH1D ("M_lljet", "M_lljet", 1000, 0, 20000);
TH1D *Muon1_pt_shape = new TH1D ("Muon1_pt_shape", "Muon1_pt_shape", 200, 0, 5000);
TH1D *Muon2_pt_shape = new TH1D ("Muon2_pt_shape", "Muon2_pt_shape", 200, 0, 5000);
TH1D *M_ll_shape = new TH1D ("M_ll_shape", "M_ll_shape", 500, 0, 8000);
TH1D *Muon1_pt = new TH1D ("Muon1_pt", "Muon1_pt", 200, 0, 5000);
TH1D *Muon2_pt = new TH1D ("Muon2_pt", "Muon2_pt", 200, 0, 5000);
TH1D *Muon1_phi = new TH1D ("Muon1_phi", "Muon1_phi", 200, -3, 3);
TH1D *Muon2_phi = new TH1D ("Muon2_phi", "Muon2_phi", 200, -3, 3);
TH1D *Muon1_eta = new TH1D ("Muon1_eta", "Muon1_eta", 200, -4, 4);
TH1D *Muon2_eta = new TH1D ("Muon2_eta", "Muon2_eta", 200, -4, 4);
TH1D *M_ll = new TH1D ("M_ll", "M_ll", 500, 0, 8000);
TH1D *BoostJet_pt = new TH1D ("BoostJet_pt", "BoostJet_pt", 200, 0, 5000);
TH1D *BoostJet_phi = new TH1D ("BoostJet_phi", "BoostJet_phi", 200, -3, 3);
TH1D *BoostJet_eta = new TH1D ("BoostJet_eta", "BoostJet_eta", 200, -4, 4);

cout << a_->GetEntries() << endl;
int tot = 0;
double deltaR = 0, deltaEta = 0, deltaPhi = 0;
double wg = 0;
double lumi = 41529;
for (Int_t i=0;i<a_->GetEntries();i++){ 
 a_->GetEntry(i);
 tot = a_->GetEntries();
 if (i%100000 == 0)cout << i << " eventi analizzati su " << tot << endl;

 if (Jet_pt->size() > 1) {
  deltaEta = fabs(Jet_eta->at(0) - Jet_eta->at(1));
  deltaPhi = Jet_phi->at(0) - Jet_phi->at(1);
  if(deltaPhi > M_PI) deltaPhi -= 2*M_PI;
  else if(deltaPhi <= -M_PI) deltaPhi +=2*M_PI;
  deltaR = sqrt(pow(deltaEta,2)+pow(deltaPhi,2));
  Jet_separation->Fill(deltaR);
  Jet_separation_vs_pt->Fill(deltaR,Jet_pt->at(0));
 }

 if (Muon_pt->size() > 1 && BoostedJet_pt->size() > 0){

  if(HLT_Mu==1 && Muon_pt->at(0) > 70)

  if (numOfHighptMu==2 && Muon_pt->at(0) > 250 && Muon_pt->at(1) > 150 && fabs(Muon_eta->at(0))<2.4 && fabs(Muon_eta->at(1))<2.4 
    && M_leplep>700 && BoostedJet_pt->at(0)>300 && numOfVetoEle==0 && numOfBoostedJets>= 1 && HLT_Mu==1){

   wg = lumi*lumi_wgt*lepsf_evt;

   //printf("wg:%f lepsf_evt:%f lumi_wgt:%f lumi:%f \n", wg, lepsf_evt, lumi_wgt, lumi);

   num->Fill(1);
   M_lljet->Fill(M_leplepBjet, wg);
   Muon1_pt_shape->Fill(Muon_pt->at(0));
   Muon2_pt_shape->Fill(Muon_pt->at(1)); 
   M_ll_shape->Fill(M_leplep,wg);
   Muon1_pt->Fill(Muon_pt->at(0),wg);
   Muon1_eta->Fill(Muon_eta->at(0),wg);
   Muon1_phi->Fill(Muon_phi->at(0),wg);
   Muon2_pt->Fill(Muon_pt->at(1),wg);
   Muon2_eta->Fill(Muon_eta->at(1),wg);
   Muon2_phi->Fill(Muon_phi->at(1),wg);
   M_ll->Fill(M_leplep,wg);
   BoostJet_pt->Fill(BoostedJet_pt->at(0),wg);
   BoostJet_eta->Fill(BoostedJet_eta->at(0),wg);
   BoostJet_phi->Fill(BoostedJet_phi->at(0),wg);
  
 
  }
 }
}

TFile *f = new TFile("HM_plot/HighMass_Selection_Sign_mumu_L11_M7000.root", "RECREATE");
num->Write();
M_lljet->Write();
Jet_separation->Write();
Jet_separation_vs_pt->Write();
Muon1_pt_shape->Write();
Muon2_pt_shape->Write();
M_ll_shape->Write();
Muon1_pt->Write();
Muon2_pt->Write();
Muon1_eta->Write();
Muon2_eta->Write();
Muon1_phi->Write();
Muon2_phi->Write();
M_ll->Write();
BoostJet_pt->Write();
BoostJet_eta->Write();
BoostJet_phi->Write();

f->Write();
f->Close();

}
 
