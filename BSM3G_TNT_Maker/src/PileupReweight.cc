#include "BSMFramework/BSM3G_TNT_Maker/interface/PileupReweight.h"
PileupReweight::PileupReweight(std::string name, TTree* tree, bool debug, const pset& iConfig, edm::ConsumesCollector && ic):
  baseTree(name,tree,debug)
{
  if(debug) std::cout<<"in PileupReweight constructor"<<std::endl;
  _MiniAODv2 = iConfig.getParameter<bool>("MiniAODv2");
  _is_data   = iConfig.getParameter<bool>("is_data");
  PUInfo_         = ic.consumes<std::vector< PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  PUReweightfile_ = iConfig.getParameter<edm::FileInPath>("MyDataPileupHistogram2018");
  //MyPUReweightfile_ = iConfig.getParameter<edm::FileInPath>("MyPUReweightfile");
  MinBiasUpReweightfile_ = iConfig.getParameter<edm::FileInPath>("MyDataPileupHistogramUp2018");
  MinBiasDownReweightfile_ = iConfig.getParameter<edm::FileInPath>("MyDataPileupHistogramDown2018");
  // Get data distribution from file
  const char *filePath = PUReweightfile_.fullPath().c_str();
  TFile file(filePath, "READ");
  TH1* h = NULL;
  file.GetObject("pileup",h);
  if( h == NULL ) {
    std::cerr << "\n\nERROR in PUWeight: Histogram 'pileup' does not exist in file \n.";
    throw std::exception();
  }
  h->SetDirectory(0);
  file.Close();

  // Get Mydata distribution from file
  /*
  const char *MyfilePath = MyPUReweightfile_.fullPath().c_str();
  TFile Myfile(MyfilePath, "READ");
  TH1* Myhist = NULL;
  Myfile.GetObject("pileup",Myhist);
  if( Myhist == NULL ) {
    std::cerr << "\n\nERROR in PUWeight: Histogram 'pileup' does not exist in file \n.";
    throw std::exception();
  }
  Myhist->SetDirectory(0);
  Myfile.Close();
  */
  // Get MydataUp distribution from file
  const char *MinBiasUpfilePath = MinBiasUpReweightfile_.fullPath().c_str();
  TFile MinBiasUpfile(MinBiasUpfilePath, "READ");
  TH1* MinBiasUphist = NULL;
  MinBiasUpfile.GetObject("pileup",MinBiasUphist);
  if( MinBiasUphist == NULL ) {
    std::cerr << "\n\nERROR in PUWeight: Histogram 'pileup' does not exist in file \n.";
    throw std::exception();
  }
  MinBiasUphist->SetDirectory(0);
  MinBiasUpfile.Close();
  
  // Get MydataDown distribution from file
  const char *MinBiasDownfilePath = MinBiasDownReweightfile_.fullPath().c_str();
  TFile MinBiasDownfile(MinBiasDownfilePath, "READ");
  TH1* MinBiasDownhist = NULL;
  MinBiasDownfile.GetObject("pileup",MinBiasDownhist);
  if( MinBiasDownhist == NULL ) {
    std::cerr << "\n\nERROR in PUWeight: Histogram 'pileup' does not exist in file \n.";
    throw std::exception();
  }
  MinBiasDownhist->SetDirectory(0);
  MinBiasDownfile.Close();
  
  // Computing weights
  // Store probabilites for each pu bin
  unsigned int nPUMax = 0;
  double *npuProbs = 0;
  
 //------------ 2018_25ns_JuneProjectionFull18 ------------//
   nPUMax =99;
   double npu_JuneProjectionFull18[nPUMax] = {4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05, 3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473, 0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138, 0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411, 0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554, 0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895, 0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877, 0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612, 0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551, 0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934, 0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915, 0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932, 0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885, 0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012, 0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05, 2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06, 3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07, 5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07, 1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08, 6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08};
  npuProbs = npu_JuneProjectionFull18;
 
  //------------------------------------------//

  // Check that binning of data-profile matches MC scenario
  if( nPUMax != static_cast<unsigned int>(h->GetNbinsX()) ) {
    std::cerr << "\n\nERROR number of bins (" << h->GetNbinsX() << ") in data PU-profile does not match number of bins (" << nPUMax << ") in MC scenario " << std::endl;
    throw std::exception();
  }
  /*
  if( nPUMax != static_cast<unsigned int>(Myhist->GetNbinsX()) ) {
    std::cerr << "\n\nERROR number of bins (" << Myhist->GetNbinsX() << ") in Mydata PU-profile does not match number of bins (" << nPUMax << ") in MC scenario " << std::endl;
    throw std::exception();
  }
  */

  if( nPUMax != static_cast<unsigned int>(MinBiasUphist->GetNbinsX()) ) {
    std::cerr << "\n\nERROR number of bins (" << MinBiasUphist->GetNbinsX() << ") in MinBiasUpdata PU-profile does not match number of bins (" << nPUMax << ") in MC scenario " << std::endl;
    throw std::exception();
  }

  if( nPUMax != static_cast<unsigned int>(MinBiasDownhist->GetNbinsX()) ) {
    std::cerr << "\n\nERROR number of bins (" << MinBiasDownhist->GetNbinsX() << ") in MinBiasDowndata PU-profile does not match number of bins (" << nPUMax << ") in MC scenario " << std::endl;
    throw std::exception();
  }

  std::vector<double> result(nPUMax,0.);
  //std::vector<double> Myresult(nPUMax,0.);
  std::vector<double> MinBiasUpresult(nPUMax,0.);
  std::vector<double> MinBiasDownresult(nPUMax,0.);
  double s = 0.;
  //double Mys = 0.;
  double MinBiasUps = 0.;
  double MinBiasDowns = 0.;
  for(unsigned int npu = 0; npu < nPUMax; ++npu) {
    const double npu_estimated = h->GetBinContent(h->GetXaxis()->FindBin(npu));
    //const double Mynpu_estimated = Myhist->GetBinContent(h->GetXaxis()->FindBin(npu));
    const double MinBiasUpnpu_estimated = MinBiasUphist->GetBinContent(h->GetXaxis()->FindBin(npu));
    const double MinBiasDownnpu_estimated = MinBiasDownhist->GetBinContent(h->GetXaxis()->FindBin(npu));
    result[npu] = npu_estimated / npuProbs[npu];
    //Myresult[npu] = Mynpu_estimated / npuProbs[npu];
    MinBiasUpresult[npu] = MinBiasUpnpu_estimated / npuProbs[npu];
    MinBiasDownresult[npu] = MinBiasDownnpu_estimated / npuProbs[npu];
    s += npu_estimated;
    //Mys += Mynpu_estimated;
    MinBiasUps += MinBiasUpnpu_estimated;
    MinBiasDowns += MinBiasDownnpu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(unsigned int npu = 0; npu < nPUMax; ++npu) {
    result[npu] /= s;
    //Myresult[npu] /= Mys;
    MinBiasUpresult[npu] /= MinBiasUps;
    MinBiasDownresult[npu] /= MinBiasDowns;
  }

  puWeigths_ = result;
 // MypuWeigths_ = Myresult;
  MinBiasUpWeigths_ = MinBiasUpresult;
  MinBiasDownWeigths_ = MinBiasDownresult;
  nPUMax_ = puWeigths_.size();
  //MynPUMax_ = MypuWeigths_.size();
  MinBiasUpnPUMax_ = MinBiasUpWeigths_.size();
  MinBiasDownnPUMax_ = MinBiasDownWeigths_.size();

  // Clean up
  delete h;
  //delete Myhist;
  delete MinBiasUphist;
  delete MinBiasDownhist;

  SetBranches();
}
PileupReweight::~PileupReweight(){
  delete tree_;
}

void PileupReweight::Fill(const edm::Event& iEvent){
  if(debug_) std::cout<<"getting PileupReweight info"<<std::endl;
  double w = 1.;
  //double Myw = 1.;
  double MinBiasUpw = 1.;
  double MinBiasDownw = 1.;
  if(!_is_data) {
    Handle<std::vector< PileupSummaryInfo > >  PUInfo;
    iEvent.getByToken(PUInfo_, PUInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    float nPU = -1;
    for(PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) { 
	nPU = PVI->getTrueNumInteractions();
	continue;
      }
    }
    if(nPU<0) nPU=0;
    if( nPU >= nPUMax_ ) {
      //std::cerr << "WARNING: Number of PU vertices = " << nPU << " out of histogram binning." << std::endl;
      // In case nPU is out-of data-profile binning,
      // use weight from last bin
      w = puWeigths_.back();
    } else {
      w = puWeigths_.at(nPU);
    }
    /*
    if( nPU >= MynPUMax_ ) {
      Myw = MypuWeigths_.back();
    } else {
      Myw = MypuWeigths_.at(nPU);
    }
    */
    if( nPU >= MinBiasUpnPUMax_ ) {
      MinBiasUpw = MinBiasUpWeigths_.back();
    } else {
      MinBiasUpw = MinBiasUpWeigths_.at(nPU);
    }
    if( nPU >= MinBiasDownnPUMax_ ) {
      MinBiasDownw = MinBiasDownWeigths_.back();
    } else {
      MinBiasDownw = MinBiasDownWeigths_.at(nPU);
    }
  }
  PUWeight=w;
 // MyPUWeight=Myw;
  MinBiasUpWeight=MinBiasUpw;
  MinBiasDownWeight=MinBiasDownw;
  if(debug_) std::cout<<"got PileupReweight info"<<std::endl;
}
void PileupReweight::SetBranches(){
  if(debug_) std::cout<<"setting branches: calling AddBranch of PileupReweight"<<std::endl;
  AddBranch(&PUWeight,"PUWeight");
  //AddBranch(&MyPUWeight,"MyPUWeight");
  AddBranch(&MinBiasUpWeight,"MinBiasUpWeight");
  AddBranch(&MinBiasDownWeight,"MinBiasDownWeight");
}
