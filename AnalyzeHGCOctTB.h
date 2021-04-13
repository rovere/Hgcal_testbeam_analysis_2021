#ifndef AnalyzeHGCOctTB_H
#define AnalyzeHGCOctTB_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

class AnalyzeHGCOctTB : public HGCNtupleVariables{

 public:
  AnalyzeHGCOctTB(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *config="alpha", const char* energy = "-1");
  ~AnalyzeHGCOctTB();
  Bool_t   FillChain(TChain *chain, TChain *chain2, TChain *chain3, const TString &inputFileList);
  /* Bool_t   FillChain(TChain *chain, const TString &inputFileList); */
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *);
  void     BookHistogram(const char *, const char *, const char* energy);

  void moduleMap_init(const char *);
  void Alignment_Map_Init();
  void Noise_Map_Init();
  void layerPosition_init();
  std::vector<bool> *noise_flag;
  TFile *oFile;
  const char *conf_;  
  int inEnergy_;
  int count=0, count_afterCuts =0;//,count_nflag=0,count_nchannel=0,count_nchanLayer1=0, count_nflag0=0;
  int count_FH=0;//, count_FH_adiMask=0,count_adiMask=0;
  int count_FH_afterCuts=0;

  // float rechits_layer[80]={};
  // float rechits_hgcallayer[40]={};
  // float	rechits_ahcallayer[39]={};
  float ratio=0;
  float Esum_rechits_EE_nocut=0;
  float Esum_rechits_EE=0;
  float Esum_rechits_EE_allcut=0;
  float	Esum_rechits_FH=0;
  float Esum_rechits_AH=0;
  float  Esum_rechits_AH_nocut=0;
  float total_rechits_energy=0;
  float total_rechits_energy_nocut=0;
  float Esum_rechits_FH_nocut=0;
  float Esum_rechits_FH_allcut=0;
  ///                           /////
  /// variables in GeV ////////
  //ratio FH = FH/(total)
  //ratio EE = EE/(total).LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL,
  float Energy_EE_inGeV,Energy_FH_inGeV,Energy_AH_inGeV,totalEnergy_inGeV,ratio_inGeV_EE,ratio_inGeV_FH,ratio_EE_inMips,ratio_FH_inMips,ratio_AH_inMips,ratio_inGeV_AH;
  /// Declare histograms here /////
  ////histograms for energy in GeV /////
  TH2F* h_Energy_EEvsFH_inGeV;
  TH2F* h_Energy_EEvsAH_inGeV;
  TH2F* h_Energy_FHvsAH_inGeV;
  TH1F* h_Energy_EE_inGeV;
  TH1F* h_Energy_FH_inGeV;
  TH1F* h_Energy_AH_inGeV;
  TH1F* h_ratioInGeV_EE;
  TH1F* h_ratioInGeV_FH;
  TH1F* h_ratioInMips_EE;
  TH1F* h_ratioInMips_FH;
  TH1F* h_ratioInGeV_AH;
  TH1F* h_ratioInMips_AH;
  TH2F* h_SEnergy_rechits_vs_layer;
  TH2F* h_EE_vs_FH;
  TH2F* h_EE_vs_AH;
  TH2F* h_FH_vs_AH;
  //2d histograms
  TH1F* h_ratioEEvsFH;
  TH2F* h_EE_vs_FH_nocut;
  TH2F* h_EE_vs_FH_allcut;
  TH2F* h_EE_vs_AH_nocut;
  TH2F* h_FH_vs_AH_nocut;
  TH1F* h_hgcal_rechits;
  TH1F* h_beamenergy;
  TH1F* h_particle;
  TH1F* h_ahcal_rechits;
  TH1F* h_SEnergy_Rechits_EE_nocut;
  TH1F* h_SEnergy_Rechits_EE_allcut;
  TH1F* h_SEnergy_Rechits_EE;  
  TH1F* h_SEnergy_Rechits_FH;
  TH1F* h_SEnergy_Rechits_FH_nocut;
  TH1F* h_SEnergy_Rechits_FH_allcut;
  TH1F* h_SEnergy_Rechits_AH;
  TH1F* h_total_rechitsEenrgy;
  TH1F* h_total_rechitsEenrgy_nocut;
  TH1F* h_totalrechits_nocut;
  TH1F* h_totalrechits_allcut;                                                    
  TH1F* h_totalrechits_FH_nocut;
  TH1F* h_totalrechits_FH_allcut;
  TH1F* h_totalEnergy_inGeV;
  //
  //   categorisng pion events
  float EnergySum_SSinEE=0,EnergySum_SSinEE_rmAHCAL=0,EnergySum_SSinFH=0,EnergySum_SSinFH_rmEE=0,EnergySum_MipLike=0,EnergySum_rmMipLike;
  float EnergySum_SSinEE_inMips=0,EnergySum_SSinEE_rmAHCAL_inMips=0,EnergySum_SSinFH_inMips=0,EnergySum_SSinFH_rmEE_inMips=0,EnergySum_MipLike_inMips=0,EnergySum_rmMipLike_inMips;
  TH1F* h_EnergySum_SSinEE;
  TH1F* h_EnergySum_SSinEE_rmAHCAL;
  TH1F* h_EnergySum_SSinFH;
  TH1F* h_EnergySum_MipLike;
  TH1F* h_EnergySum_SSinFH_rmEE;
  TH1F* h_EnergySum_rmMipLike;
   TH1F* h_EnergySum_SSinEE_inMips;
  TH1F* h_EnergySum_SSinEE_rmAHCAL_inMips;
  TH1F* h_EnergySum_SSinFH_inMips;
  TH1F* h_EnergySum_MipLike_inMips;
  TH1F* h_EnergySum_SSinFH_rmEE_inMips;
  TH1F* h_EnergySum_rmMipLike_inMips;

};
#endif

#ifdef AnalyzeHGCOctTB_cxx

void AnalyzeHGCOctTB::BookHistogram(const char *outFileName, const char* conf,  const char* energy) {

  conf_ = conf;
  oFile = new TFile(outFileName, "recreate");

  int en = atoi(energy);
  float xlow = 0.0;
  float xhigh = en*100*1.5;
  //int xbin = (xhigh-xlow)/6;
  int xbin = 500;//*0.5*en;
  ////////// Book histograms here  //////////////
  //h_rechit_layer= new TH1F("h_rechit_layer","h_rechit_layer",80,0,80);replace string
  h_EnergySum_SSinEE=new TH1F("h_EnergySum_SSinEE","h_EnergySum_SSinEE",xbin,xlow,3.0*en);
  h_EnergySum_SSinEE_rmAHCAL=new TH1F("h_EnergySum_SSinEE_rmAHCAL","h_EnergySum_SSinEE_rmAHCAL",xbin,xlow,3.0*en);
  h_EnergySum_SSinFH=new TH1F("h_EnergySum_SSinFH","h_EnergySum_SSinFH",xbin,xlow,3.0*en);
  h_EnergySum_SSinFH_rmEE=new TH1F("h_EnergySum_SSinFH_rmEE","h_EnergySum_SSinFH_rmEE",xbin,xlow,3.0*en);
  h_EnergySum_MipLike=new TH1F("h_EnergySum_MipLike","h_EnergySum_MipLike",xbin,xlow,3.0*en);
  h_EnergySum_rmMipLike=new TH1F("h_EnergySum_rmMipLike","h_EnergySum_rmMipLike",xbin,xlow,3.0*en);


   h_EnergySum_SSinEE_inMips=new TH1F("h_EnergySum_SSinEE_inMips","h_EnergySum_SSinEE_inMips",xbin,xlow,xhigh);
  h_EnergySum_SSinEE_rmAHCAL_inMips=new TH1F("h_EnergySum_SSinEE_rmAHCAL_inMips","h_EnergySum_SSinEE_rmAHCAL_inMips",xbin,xlow,xhigh);
  h_EnergySum_SSinFH_inMips=new TH1F("h_EnergySum_SSinFH_inMips","h_EnergySum_SSinFH_inMips",xbin,xlow,xhigh);
  h_EnergySum_SSinFH_rmEE_inMips=new TH1F("h_EnergySum_SSinFH_rmEE_inMips","h_EnergySum_SSinFH_rmEE_inMips",xbin,xlow,xhigh);
  h_EnergySum_MipLike_inMips=new TH1F("h_EnergySum_MipLike_inMips","h_EnergySum_MipLike_inMips",xbin,xlow,xhigh);
  h_EnergySum_rmMipLike_inMips=new TH1F("h_EnergySum_rmMipLike_inMips","h_EnergySum_rmMipLike_inMips",xbin,xlow,xhigh);

  h_Energy_EEvsFH_inGeV=new TH2F("h_Energy_EEvsFH_inGeV","h_Energy_EEvsFH_inGeV",xbin,xlow,3.0*en,xbin,xlow,3.0*en);
  h_Energy_EEvsAH_inGeV=new TH2F("h_Energy_EEvsAH_inGeV","h_Energy_EEvsAH_inGeV",xbin,xlow,3.0*en,xbin,xlow,3.0*en);
  h_Energy_FHvsAH_inGeV=new TH2F("h_Energy_FHvsAH_inGeV","h_Energy_FHvsAH_inGeV",xbin,xlow,3.0*en,xbin,xlow,3.0*en);
  h_Energy_EE_inGeV=new TH1F("h_Energy_EE_inGeV","h_Energy_EE_inGeV",xbin,xlow,3.0*en);
  h_Energy_FH_inGeV=new TH1F("h_Energy_FH_inGeV","h_Energy_FH_inGeV",xbin,xlow,3.0*en);
  h_Energy_AH_inGeV=new TH1F("h_Energy_AH_inGeV","h_Energy_AH_inGeV",xbin,xlow,3.0*en);

  
  h_ratioEEvsFH= new TH1F("h_ratioEEvsFH","h_ratioEEvsFH",xbin,xlow,3*xhigh);
  h_ratioInGeV_EE=new TH1F("h_ratioInGeV_EE","EE/(FH+EE+AH)",xbin,xlow,2.0);
  h_ratioInGeV_FH=new TH1F("h_ratioInGeV_FH","FH/(FH+EE+AH)",xbin,xlow,2.0);
  h_ratioInGeV_AH=new TH1F("h_ratioInGeV_AH","AH/(FH+EE+AH)",xbin,xlow,2.0);
  h_ratioInMips_EE=new TH1F("h_ratioInMips_EE","EE/(FH+EE+AH)",xbin,xlow,2.0);
  h_ratioInMips_FH=new TH1F("h_ratioInMips_FH","FH/(FH+EE+AH)",xbin,xlow,2.0);
  h_ratioInMips_AH=new TH1F("h_ratioInMips_AH","AH/(FH+EE+AH)",xbin,xlow,2.0);
  h_totalEnergy_inGeV=new TH1F("h_totalEnergy_inGeV","h_totalEnergy_inGeV",xbin,xlow,3.0*en);
  h_SEnergy_rechits_vs_layer = new TH2F("h_SEnergy_rechits_vs_layer","h_rechits_vs_layer",40,0,40,500,0,5000);
  
  h_EE_vs_AH=new TH2F("h_EE_vs_AH","total rechits energy in EE vs AH",xbin,xlow,xhigh,xbin,xlow,xhigh);
  
  h_EE_vs_FH=new TH2F("h_EE_vs_FH","total rechits energy in EE vs FH",xbin,xlow,xhigh,xbin,xlow,xhigh);
  h_EE_vs_FH_nocut=new TH2F("h_EE_vs_FH_nocut","h_EE_vs_FH_nocut",xbin,xlow,xhigh,xbin,xlow,xhigh);
  h_EE_vs_FH_allcut=new TH2F("h_EE_vs_FH_allcut","h_EE_vs_FH_allcut",xbin,xlow,xhigh,xbin,xlow,xhigh);
  h_EE_vs_AH_nocut=new TH2F("h_EE_vs_AH_nocut","total rechits energy in EE vs AH no cut",xbin,xlow,xhigh,xbin,xlow,xhigh);
  h_FH_vs_AH_nocut=new TH2F("h_FH_vs_AH_nocut","total rechits energy in FH vs AH no cut",xbin,xlow,xhigh,xbin,xlow,xhigh);
  h_FH_vs_AH=new TH2F("h_FH_vs_AH","total rechits energy in FH vs AH",xbin,xlow,xhigh,xbin,xlow,xhigh);
  
  h_SEnergy_Rechits_AH=new TH1F("h_SEnergy_Rechits_AH","h_SEnergy_Rechits_AH",xbin,xlow,xhigh);
  //h_SEnergy_Rechits_AH=new TH1F("h_SEnergy_Rechits_AH","h_SEnergy_Rechits_AH",xbin,xlow,xhigh);

  h_SEnergy_Rechits_FH=new TH1F("h_SEnergy_Rechits_FH","h_SEnergy_Rechits_FH",xbin,xlow,xhigh);
  h_SEnergy_Rechits_EE_nocut=new TH1F("h_SEnergy_Rechits_EE_nocut","h_SEnergy_Rechits_EE_nocut",xbin,xlow,xhigh);
  h_SEnergy_Rechits_EE_allcut=new TH1F("h_SEnergy_Rechits_EE_allcut","h_SEnergy_Rechits_EE_allcut",xbin,xlow,xhigh);
  h_SEnergy_Rechits_EE=new TH1F("h_SEnergy_Rechits_EE","h_SEnergy_Rechits_EE",xbin,xlow,xhigh);
  
  /////////////////////
  h_SEnergy_Rechits_FH_nocut=new TH1F("h_SEnergy_Rechits_FH_nocut","h_SEnergy_Rechits_FH_nocut",xbin,xlow,xhigh);
  h_SEnergy_Rechits_FH_allcut=new TH1F("h_SEnergy_Rechits_FH_allcut","h_SEnergy_Rechits_FH_allcut",xbin,xlow,xhigh);
  ///////////////////
  
  h_total_rechitsEenrgy=new TH1F("h_total_rechitsEenrgy","h_total_rechitsEenrgy",xbin,xlow,xhigh);
  h_total_rechitsEenrgy_nocut=new TH1F("h_total_rechitsEenrgy_nocut","h_total_rechitsEenrgy_nocut",xbin,xlow,xhigh);

  h_hgcal_rechits = new TH1F("h_hgcal_rechits","h_hgcal_rechits",1500,0,1500);
  h_ahcal_rechits = new TH1F("h_ahcal_rechits","h_ahcal_rechits",1500,0,1500);
  h_beamenergy = new TH1F("h_beamenergy","h_beamenergy",320,0,320);
  h_particle = new TH1F("h_particle","h_particle",600,-300,300);
  h_totalrechits_nocut= new TH1F("h_totalrechits_nocut","h_totalrechits_nocut",1000,0,5000);
  h_totalrechits_allcut=new TH1F("h_totalrechits_allcut","h_totalrechits_allcut",1000,0,5000);
  //////
   h_totalrechits_FH_nocut= new TH1F("h_totalrechits_FH_nocut","h_totalrechits_FH_nocut",1000,0,5000);
  h_totalrechits_FH_allcut=new TH1F("h_totalrechits_FH_allcut","h_totalrechits_FH_allcut",1000,0,5000);
  
}

void AnalyzeHGCOctTB::Alignment_Map_Init() {
  /// Alignment map for config 1, it needs to be changed accorfing to configurations
  
  char* f_name = new char[200];
  sprintf(f_name,"../Alignment_Map.txt");
  std::ifstream in(f_name);
  //int layer;
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    //return;
    exit(0);
  }

  std::pair<float,float> dx_dy;
  std::pair<int, std::pair<float,float> > temp;
  int layer;
  float dx,dy;
  while(in>>layer>>dx>>dy) {
    //run_layer = std::make_pair(run,layer);
    dx_dy = std::make_pair(dx,dy);
    temp = std::make_pair(layer,dx_dy);
    align_map.insert(temp);
  }

  std::cout<<"INFO: Alignment MAP initialized successfully!!!"<<endl;
}

void AnalyzeHGCOctTB::Noise_Map_Init() {
  //Noise map for config - 1 , needs to be changed accordingly
  
  char* f_name = new char[200];
  sprintf(f_name,"../Noise_Map.txt");
  std::ifstream in(f_name);
  
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    exit(0);
  }
  std::pair<int,int> mod_chip;
  std::pair<std::pair<int,int>, float> temp;
  int layer,mod_id,mod_pos,chip;
  float noise;
  while(in>>layer>>mod_id>>mod_pos>>chip>>noise) {
    mod_chip = std::make_pair(mod_id,chip);
    temp = std::make_pair(mod_chip,noise);
    noise_map.insert(temp);
  }

  std::cout<<"INFO: Noise MAP initialized successfully!!!"<<endl;
}

void AnalyzeHGCOctTB::moduleMap_init(const char* config) {
  char *f_name = new char[200];

  if(strcmp(config,"alpha")==0 || strcmp(config,"config1")==0) {
    sprintf(f_name,"../config_maps/moduleMAP_config1.txt");
    cout<<"\n\nINFO: Mapping module configuration ALPHA (oct10-oct17) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else if(strcmp(config,"bravo")==0 || strcmp(config,"config2")==0) {
    sprintf(f_name,"../config_maps/moduleMAP_config2.txt");
    cout<<"\n\nINFO: Mapping module configuration BRAVO (17oct-22oct) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else if(strcmp(config,"charlie")==0  || strcmp(config,"config3")==0) {
    sprintf(f_name,"../config_maps/moduleMAP_config3.txt");
    cout<<"\n\nINFO: Mapping module configuration CHARLIE (23Oct-4Nov) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else {
    cout<<"\n\nERROR: Incorrect configuration entered "<<endl;
    cout<<" Allowed configuration :\n alpha = Configuration 1 (10Oct-17Nov) \n bravo = Configuration 2 (17Oct-22Oct) \n charlie = Configuration 3 (23Oct-17Nov)"<<endl;
    return;
    
  }

  std::ifstream in(f_name);
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    //return;
    exit(0);
  }


  int modID_, part_, layer_, pos_;
  cout<<"File name = "<<f_name<<endl;
  while(in>>modID_>>part_>>layer_>>pos_){
    std::pair<int, std::vector<int>> temp_pair;
    std::vector<int> temp_vector;
    temp_vector.push_back(part_);
    temp_vector.push_back(layer_);
    temp_vector.push_back(pos_);
    temp_pair = std::make_pair(modID_,temp_vector);
    module_map.insert(temp_pair);
  }

  cout<<"INFO: Module Mapping Done!!! "<<endl<<endl;


}


void AnalyzeHGCOctTB::layerPosition_init() {
  // Layer postions for config - 1, needs to be changed accordingly for other configs
  
  char *f_name = new char[200];
  sprintf(f_name,"../config1_lengths.txt");

  std::ifstream in(f_name);
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    exit(0);
  }
  int layer_;
  float  cm_, x0, nuc_int_, pi_int_;

  cout<<"File name = "<<f_name<<endl;
  while(in>>layer_>>cm_>>x0>>nuc_int_>>pi_int_){
    std::pair<int, std::vector<float>> temp_pair;
    std::vector<float> temp_vector;
    temp_vector.push_back(cm_/10.0);
    temp_vector.push_back(x0);
    temp_vector.push_back(nuc_int_);
    temp_vector.push_back(pi_int_);

    temp_pair = std::make_pair(layer_,temp_vector);
    layer_positions.insert(temp_pair);
  }

  cout<<"INFO: Layer Position Mapping Done!!! "<<endl<<endl;


}


AnalyzeHGCOctTB::AnalyzeHGCOctTB(const TString &inputFileList, const char *outFileName, const char* dataset, const char* config, const char* energy) {

  TChain *tree = new TChain("rechitntupler/hits");
  TChain *tree2 = new TChain("trackimpactntupler/impactPoints");
  TChain *tree3 = new TChain("bigtree");


  if( ! FillChain(tree, tree2, tree3, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }


  HGCNtupleVariables::Init(tree, tree2, tree3);

  BookHistogram(outFileName, config, energy);
  moduleMap_init(config);
  // Alignment_Map_Init();
  // Noise_Map_Init();
  // layerPosition_init();
  
}
Bool_t AnalyzeHGCOctTB::FillChain(TChain *chain, TChain *chain2, TChain *chain3, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
    chain2->Add(buffer.c_str());
    chain3->Add(buffer.c_str());

  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  std::cout << "No. of Entries in chain3 : " << chain3->GetEntries() << std::endl;

  return kTRUE;
}

Long64_t AnalyzeHGCOctTB::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  if (!fChain2) return -5;
  Long64_t centry2 = fChain2->LoadTree(entry);
  if (centry2 < 0) return centry2;
  if (!fChain2->InheritsFrom(TChain::Class()))  return centry2;
  TChain *chain2 = (TChain*)fChain2;
  if (chain2->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  if (!fChain3) return -5;
  Long64_t centry3 = fChain3->LoadTree(entry);
  if (centry3 < 0) return centry3;
  if (!fChain3->InheritsFrom(TChain::Class()))  return centry3;
  TChain *chain3 = (TChain*)fChain3;
  if (chain3->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  
  //if (centry==centry2)
  return centry;
  // cout<<"centry = "<<centry<<endl;
  // if(centry>0)
  //   return centry;
  // else return -1;
}

AnalyzeHGCOctTB::~AnalyzeHGCOctTB() { 

  // if (!fChain || !fChain2) return;
  // delete fChain->GetCurrentFile();
  // delete fChain2->GetCurrentFile();
  // oFile->cd();
  // oFile->Write();
  // oFile->Close();


  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif

/*  LocalWords:  Nrechit EE R1 FH GetXaxis SetTitle Sumw2 TH2F reg3 NRechits
 */
/*  LocalWords:  GetYaxis SetTitleOffset
 */
