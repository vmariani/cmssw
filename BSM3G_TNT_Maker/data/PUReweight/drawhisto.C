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

TFile *f00 = new TFile("MC_JuneProjectionFull18_configuration.root");
TFile *f01 = new TFile("MyDataPileupHistogramDown2018.root");

double bin = 1;
TString name = "pileup";
TH1F *da= (TH1F*) f00->Get(name);
TH1F *da1= (TH1F*) f01->Get(name);

TH1F *d= (TH1F*) da->Clone(); 
TH1F *d1= (TH1F*) da1->Clone();

d->SetLineColor(1); //1
d->SetLineWidth(3);
d->SetLineStyle(1);
   
d1->SetLineColor(2); //2
d1->SetLineWidth(2);
d1->SetLineStyle(1);	

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
d->GetXaxis()->SetTitle("PU");
d->GetYaxis()->SetLabelSize(0.045);
d->GetYaxis()->SetTitleSize(0.06);
double bincontent = d->GetBinContent(0);

d->Sumw2();
double w = d->Integral();//GetEntries();
d->Scale(1./w);
d1->Sumw2();
double w1 = d1->Integral();//GetEntries();
d1->Scale(1./w1);
d->Rebin(bin);
d1->Rebin(bin);

const int nbin = d->GetNbinsX();
cout << nbin << endl;

double weight[100];

for(int i= 0; i < nbin; i++){
 weight[i] = 1; 
 if (d1->GetBinContent(i) != 0){
  weight[i]=(d1->GetBinContent(i)/d->GetBinContent(i));
  cout << "bin " << i << " data: " << d1->GetBinContent(i) << " mc: " << d->GetBinContent(i)  <<  " peso " << weight[i] << " bin center " << d->GetBinCenter(i) <<  endl;
 }
}

for(int i= 0; i < nbin; i++){
 cout << weight[i] << " , " ;
} 
    
//d1->Draw("Ehist");
d->Draw("hist");
d1->Draw("histsame");
//d2->Draw("Ehistsame");

TLatex latexLabel2;
latexLabel2.SetTextSize(0.04);
latexLabel2.SetTextFont(32);
latexLabel2.SetNDC();
//latexLabel2.DrawLatex(0.4, 0.93, "CMS Preliminary, pp L = 2.52 fb^{-1} at #sqrt{s} = 13 TeV");
latexLabel2.DrawLatex(0.4, 0.93, "");
//latexLabel2.DrawLatex(0.4, 0.93, "CMS Simulations, pp at #sqrt{s} = 13 TeV");//"CMS Preliminary, pp L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
   
// poi disegno la legenda  
 
TLegend *pl = new TLegend(0.6,0.75,0.95,0.85); //0.4,0.20,0.75,0.30
pl->SetTextSize(0.03); 
pl->SetFillColor(0);
pl-> SetNColumns(1);//per piu plot sovrapposti
pl->AddEntry(d, "MC", "Column 1 line 1");// "PL");
pl->AddEntry(d1, "Data", "Column 1 line 2");// "L");
//pl->AddEntry(d2, "MC MinBias", "Column 1 line 2");
//pl->AddEntry(d7, "MC", "Column 1 line 2");
     // metti la descrizione che vuoi appaia e poi con PL intendi che prende come riferimento i point e la line style dell'histo

pl->Draw();
    //gPad->RedrawAxis(); // facoltativo, solo se l'histo finisce sotto un asse
    c1->SaveAs("PU_down.png");

    // salva immagine	
//}
    
}
