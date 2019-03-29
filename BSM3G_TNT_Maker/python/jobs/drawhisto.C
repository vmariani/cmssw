{

////////// PRIMA PARTE PER CONFIG HISTO IN CMS STYLE
  gROOT->Reset();
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);
  gStyle->SetTitleX(0.5); //title X location 
  gStyle->SetTitleY(0.96); //title Y location 
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
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.23);
  tdrStyle->SetPadRightMargin(0.05);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.05);
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

  // Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();


///////////////// FINE CONFIG



  using namespace std;
TH1F * gHisto ;

TFile *f00 = new TFile("LowMass_Selection_HLT_or_Mu100_TTpST.root");
TFile *f01 = new TFile("LowMass_Selection_HLT_or_Mu100_DY.root");
TFile *f02 = new TFile("LowMass_Selection_HLT_or_Mu100_Sign_mumu_L11_M500.root");
TFile *f03 = new TFile("LowMass_Selection_HLT_or_Mu100_Sign_mumu_L20_M1000.root");
TFile *f04 = new TFile("LowMass_Selection_TTpST.root");
TFile *f05 = new TFile("LowMass_Selection_DY.root");
TFile *f06 = new TFile("LowMass_Selection_Sign_mumu_L11_M500.root");
TFile *f07 = new TFile("LowMass_Selection_Sign_mumu_L20_M1000.root");
//MyDataDleupHistogram.root

double bin = 10;
TString name = "M_ll_shape";
TH1F *da= (TH1F*) f00->Get(name);
TH1F *da1= (TH1F*) f01->Get(name);
TH1F *da2= (TH1F*) f02->Get(name);
TH1F *da3= (TH1F*) f03->Get(name);
TH1F *da4= (TH1F*) f04->Get(name);
TH1F *da5= (TH1F*) f05->Get(name);
TH1F *da6= (TH1F*) f06->Get(name);
TH1F *da7= (TH1F*) f07->Get(name);
TH1F *d= (TH1F*) da->Clone(); 
TH1F *d1= (TH1F*) da1->Clone();
TH1F *d2= (TH1F*) da2->Clone();
TH1F *d3= (TH1F*) da3->Clone();
TH1F *d4= (TH1F*) da4->Clone();
TH1F *d5= (TH1F*) da5->Clone();
TH1F *d6= (TH1F*) da6->Clone();
TH1F *d7= (TH1F*) da7->Clone();

d->SetLineColor(1);
d->SetLineWidth(3);
d->SetLineStyle(1);
   
d1->SetLineColor(2);
d1->SetLineWidth(2);
d1->SetLineStyle(1);	

d2->SetLineColor(4);
d2->SetLineWidth(2);
d2->SetLineStyle(1);

d3->SetLineColor(6);
d3->SetLineWidth(2);
d3->SetLineStyle(1);

d4->SetLineColor(8);
d4->SetLineWidth(2);
d4->SetLineStyle(1);
    
d5->SetLineColor(93);
d5->SetLineWidth(2);
d5->SetLineStyle(1);

d6->SetLineColor(46);
d6->SetLineWidth(3);
d6->SetLineStyle(1);

d7->SetLineColor(28);
d7->SetLineWidth(3);
d7->SetLineStyle(1);

TCanvas* c1 = new TCanvas("c1","c1",0,0,800,600);
TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.99);
c1_1->Draw();
c1_1->cd();
c1_1->SetTopMargin(0.1);
c1_1->SetBottomMargin(0.15);
c1_1->SetRightMargin(0.02);
c1_1->SetLeftMargin(0.15);

d->GetXaxis()->SetTitleOffset(1.);
d->GetXaxis()->SetLabelSize(0.045);
d->GetXaxis()->SetTitleSize(0.06);
d->GetYaxis()->SetTitleOffset(1.2);
d->GetYaxis()->SetTitle("Events");
//"<d^{2}N_{ch} / (d#pt d#pt)>");  sono esempi di modi di utilizzare i caratteri latex
d->GetXaxis()->SetTitle("M_{ll}");
d2->GetXaxis()->SetTitle("M_{ll}");
//d->GetXaxis()->SetTitle("m_{ll} [GeV]");
d->GetYaxis()->SetLabelSize(0.045);
d->GetYaxis()->SetTitleSize(0.06);
double bincontent = d->GetBinContent(0);

d->Sumw2();
double w = d->Integral();//GetEntries();
d->Scale(1./w);
d1->Sumw2();
double w1 = d1->Integral();//GetEntries();
d1->Scale(1./w1);
//d->SetMaximum(0.6);
d2->Sumw2();
double w2 = d2->Integral();//GetEntries();
d2->Scale(1./w2);
d3->Sumw2();
double w3 = d3->Integral();//GetEntries();
d3->Scale(1./w3);
d4->Sumw2();
double w4 = d4->Integral();//GetEntries();
d4->Scale(1./w4);
d5->Sumw2();
double w5 = d5->Integral();//GetEntries();
d5->Scale(1./w5);
d6->Sumw2();
double w6 = d6->Integral();
d6->Scale(1./w6);
d7->Sumw2();
double w7 = d7->Integral();
d7->Scale(1./w7);

d->Rebin(bin);
d1->Rebin(bin);
d2->Rebin(bin);
d3->Rebin(bin);
d4->Rebin(bin);
d5->Rebin(bin);
d6->Rebin(bin);
d7->Rebin(bin);

//d->Draw("Ehist");
//d1->Draw("Ehistsame");
d2->Draw("Ehistsame");
d3->Draw("Ehistsame");
//d4->Draw("Ehistsame");
//d5->Draw("Ehistsame");
d6->Draw("Ehistsame");
d7->Draw("Ehistsame");

TLatex latexLabel2;
latexLabel2.SetTextSize(0.04);
latexLabel2.SetTextFont(32);
latexLabel2.SetNDC();
latexLabel2.DrawLatex(0.4, 0.93, "");
   
// poi disegno la legenda  
 
TLegend *pl = new TLegend(0.48,0.8,0.98,0.95); //0.4,0.20,0.75,0.30
pl->SetTextSize(0.03); 
pl->SetFillColor(0);
pl-> SetNColumns(1);//per piu plot sovrapposti
//pl->AddEntry(d, "TT+ST HLTMu50 || pT_{#mu1} > 100 GeV", "Column 1 line 1");// "PL");
//pl->AddEntry(d1, "DY HLTMu50 || pT_{#mu1} > 100 GeV", "Column 1 line 2");// "L");
pl->AddEntry(d2, "L11 M500 HLTMu50 || pT_{#mu1} > 100 GeV", "Column 1 line 2");
pl->AddEntry(d3, "L20 M1000 HLTMu50 || pT_{#mu1} > 100 GeV", "Column 1 line 1");
//pl->AddEntry(d4, "TT+ST HLTMu50", "Column 1 line 1");
//pl->AddEntry(d5, "DY HLTMu50", "Column 1 line 1");
pl->AddEntry(d6, "L11 M500 HLTMu50", "Column 1 line 1");
pl->AddEntry(d7, "L20 M1000 HLTMu50", "Column 1 line 1");

     // metti la descrizione che vuoi appaia e poi con PL intendi che prende come riferimento i point e la line style dell'histo

pl->Draw();
    //gPad->RedrawAxis(); // facoltativo, solo se l'histo finisce sotto un asse
    c1->SaveAs("Immagini/LowMass_Selection_HLT_"+name+"_sgn.png");
    c1->SaveAs("Immagini/LowMass_Selection_HLT_"+name+"_sgn.C");

    // salva immagine	
//}
    
}
