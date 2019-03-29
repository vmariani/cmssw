{
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  //gStyle->SetTitleX(1); //title X location 
  //gStyle->SetTitleY(2); //title Y location 
  gStyle->SetPaintTextFormat(".2f");

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(true);
  tdrStyle->SetGridColor(1);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistFillColor(63);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

//  tdrStyle->SetEndErrorSize(0);
  tdrStyle->SetErrorX(0.);
//  tdrStyle->SetErrorMarker(20);
  
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(2);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.15);
  tdrStyle->SetPadBottomMargin(0.138);
  tdrStyle->SetPadLeftMargin(0.23);
  tdrStyle->SetPadRightMargin(0.05);

  // For the Global title:

  //tdrStyle->SetOptTitle(0);
  //tdrStyle->SetTitleFont(42);
  //tdrStyle->SetTitleColor(1);
  //tdrStyle->SetTitleTextColor(1);
  //tdrStyle->SetTitleFillColor(10);
  //tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  //tdrStyle->SetTitleColor(1, "XYZ");
  //tdrStyle->SetTitleFont(42, "XYZ");
  //tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  //tdrStyle->SetTitleXOffset(0.9);
  //tdrStyle->SetTitleYOffset(1.05);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  tdrStyle->cd();



  using namespace std;


vector<TString> name;   vector<TString> axis; vector<int> rebin;  vector<int> logy; vector<double> mini; vector<double> maxi;

name.push_back("Muon_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(25); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("Muon_eta"); axis.push_back("#eta"); rebin.push_back(20); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("Muon_phi"); axis.push_back("#phi"); rebin.push_back(20); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("patElectron_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(25); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("patElectron_eta"); axis.push_back("#eta"); rebin.push_back(20); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("patElectron_phi"); axis.push_back("#phi"); rebin.push_back(20); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("Jet_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(25); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("Jet_eta"); axis.push_back("#eta"); rebin.push_back(25); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("Jet_phi"); axis.push_back("#phi"); rebin.push_back(25); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("BoostedJet_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(50); logy.push_back(0); mini.push_back(100); maxi.push_back(500);

name.push_back("BoostedJet_eta"); axis.push_back("#eta"); rebin.push_back(20); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("BoostedJet_phi"); axis.push_back("#phi"); rebin.push_back(20); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("numOfHighptMu"); axis.push_back("# hight pT #mu"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfLooseMu"); axis.push_back("# loose #mu"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfMediumMu"); axis.push_back("# medium #mu"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfTightMu"); axis.push_back("# tight pT #mu"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfHighptEle"); axis.push_back("# hight pT ele"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfLooseEle"); axis.push_back("# loose ele"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfMediumEle"); axis.push_back("# medium ele"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfTightEle"); axis.push_back("# tight pT ele"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfJets"); axis.push_back("# jets"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);

name.push_back("numOfBoostedJets"); axis.push_back("# boosted jets"); rebin.push_back(10); logy.push_back(0); mini.push_back(0); maxi.push_back(10);
/*
name.push_back("Muon1_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(2); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("Muon1_eta"); axis.push_back("#eta"); rebin.push_back(4); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("Muon1_phi"); axis.push_back("#phi"); rebin.push_back(4); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("Muon2_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(2); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("Muon2_eta"); axis.push_back("#eta"); rebin.push_back(4); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("Muon2_phi"); axis.push_back("#phi"); rebin.push_back(4); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("Ele1_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(2); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("Ele1_eta"); axis.push_back("#eta"); rebin.push_back(4); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("Ele1_phi"); axis.push_back("#phi"); rebin.push_back(4); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("Ele2_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(2); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("Ele2_eta"); axis.push_back("#eta"); rebin.push_back(4); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("Ele2_phi"); axis.push_back("#phi"); rebin.push_back(4); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("BoostedJet1_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(2); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("BoostedJet1_eta"); axis.push_back("#eta"); rebin.push_back(4); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("BoostedJet1_phi"); axis.push_back("#phi"); rebin.push_back(4); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("BoostedJet2_pt"); axis.push_back("p_T [GeV]"); rebin.push_back(2); logy.push_back(0); mini.push_back(0); maxi.push_back(100);

name.push_back("BoostedJet2_eta"); axis.push_back("#eta"); rebin.push_back(4); logy.push_back(0); mini.push_back(-3); maxi.push_back(3);

name.push_back("BoostedJet2_phi"); axis.push_back("#phi"); rebin.push_back(4); logy.push_back(0); mini.push_back(-4); maxi.push_back(4);

name.push_back("M_ele1ele2Bjet1"); axis.push_back("M(e1e2bj1)"); rebin.push_back(4); logy.push_back(0);mini.push_back(0); maxi.push_back(500);

name.push_back("M_ele1ele2"); axis.push_back("M(e1e2)"); rebin.push_back(4); logy.push_back(0);mini.push_back(0); maxi.push_back(500);

name.push_back("M_elemu"); axis.push_back("M(e#mu)"); rebin.push_back(4); logy.push_back(0);mini.push_back(0); maxi.push_back(500);

name.push_back("M_elemuBjet1"); axis.push_back("M(e#mubj1)"); rebin.push_back(4); logy.push_back(0);mini.push_back(0); maxi.push_back(500);

name.push_back("M_mu1mu2"); axis.push_back("M(#mu1#mu2)"); rebin.push_back(4); logy.push_back(0);mini.push_back(0); maxi.push_back(500);

name.push_back("M_mu1mu2Bjet1"); axis.push_back("M(#mu1#mu2bj1)"); rebin.push_back(4); logy.push_back(0);mini.push_back(0); maxi.push_back(500);
*/
TFile *f00 = new TFile("plot_ZZ.root");
TFile *f01 = new TFile("plot_ZZ_2016.root");

 TTree *data_t=(TTree*) f00->Get("BOOM");
 TTree *data1_t=(TTree*) f01->Get("BOOM");

 for(int i=0; i<name.size(); i++){

 cout << name[i]+">>+data_"+name[i]+"("+std::to_string(rebin[i])+","+std::to_string(mini[i])+","+std::to_string(maxi[i])+")" << endl;
 data_t->Draw(name[i]+">>+data_"+name[i]+"("+std::to_string(rebin[i])+","+std::to_string(mini[i])+","+std::to_string(maxi[i])+")");
 TH1F *d = (TH1F*)gDirectory->Get("data_"+name[i]);

 data1_t->Draw(name[i]+">>+data1_"+name[i]+"("+std::to_string(rebin[i])+","+std::to_string(mini[i])+","+std::to_string(maxi[i])+")");
 TH1F *d1 = (TH1F*)gDirectory->Get("data1_"+name[i]);


 
 //d->Rebin(rebin[i],"",xbins);
 //tdrStyle->SetOptLogy(logy[i]);
 //d1->Rebin(rebin[i],"",xbins);
  
 d->SetLineColor(kBlue);
 d->SetFillStyle(3001);
 d->SetLineWidth(3);
 d->SetLineStyle(1); 

 d1->SetLineColor(kRed);
 d1->SetFillStyle(3001);
 d1->SetLineWidth(3);
 d1->SetLineStyle(1); 

 TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
	
 TPad *c1_1 = new TPad("c1_1", "newpad",0.0,0.01,0.99,0.99);
 c1_1->Draw();
 c1_1->cd();
 c1_1->SetTopMargin(0.1);
 c1_1->SetBottomMargin(0.15);
 c1_1->SetRightMargin(0.02);
 c1_1->SetLeftMargin(0.15);
   // nTrans_SisCon5->Draw();	

 d->GetXaxis()->SetTitleOffset(1.);
 d->GetXaxis()->SetLabelSize(0.045);
 d->GetXaxis()->SetTitleSize(0.06);
 d->GetYaxis()->SetTitleOffset(1.2);
 d->GetYaxis()->SetTitle("# events");
 d->GetXaxis()->SetTitle(axis[i]);
 d->GetYaxis()->SetLabelSize(0.045);
 d->GetYaxis()->SetTitleSize(0.06);
 d->GetXaxis()->SetRangeUser(mini[i],maxi[i]); 		
 //d->SetMinimum(0.0000001);//d->GetMinimum()+.000001);

 d->Sumw2();
 d1->Sumw2();
 double w = d->Integral();
// double w = data_t->GetEntries();
 cout << w << endl;
 d->Scale(1./w);

 double w1 = d1->Integral(); 
// double w1 = data1_t->GetEntries();
 d1->Scale(1./w1);

 d->SetMaximum(d1->GetMaximum()*2);

 d->Draw("Ehist");	
 d1->Draw("Ehistsame");

 TLatex latexLabel2;
 latexLabel2.SetTextSize(0.04);
 latexLabel2.SetTextFont(32);
 latexLabel2.SetNDC();
    //latexLabel2.DrawLatex(0.4, 0.93, "CMS Simulations");//"CMS Preliminary, pp L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
 TLegend *pl = new TLegend(0.6,0.75,0.75,0.85);
 pl->SetTextSize(0.03); 
 pl->SetFillColor(0);
 pl-> SetNColumns(1);
 pl->AddEntry(d, "2017", "Column 1 line 1");
 pl->AddEntry(d1, "2016", "Column 1 line 2");
 pl->Draw();
 gPad->RedrawAxis();
 c1->SaveAs("Immagini/ZZ_"+name[i]+"_2016_vs_2017.png");
   
 }
}
//gSystem->Exit(0); 
