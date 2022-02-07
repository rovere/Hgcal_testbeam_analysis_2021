#define AnalyzeHGCOctTB_cxx
#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>
#include<TF1.h>
using namespace std;



// chip 3022,44,3028




int main(int argc, char* argv[])//, int argvv[])
{

  if (argc < 6) {
    cerr << "Please give 7 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" " << "configuration" 
	 <<" " << "chi ethod" << "folderID"<< endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *config          = argv[4];
  const char *chi_method = argv[5];
  const char *folderID = argv[6];
  //  const char *max_ = argv[7];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, chi_method,folderID);
  cout << "dataset " << data << " " << endl;
  cout << "configuration " << config << " " << endl;
  cout << "energy " <<chi_method << " " << endl;
 int chi2_method = atoi(chi_method);
  hgcOctTB.EventLoop(data,chi_method,folderID);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data, const char *chi_method, const char *folderID) {


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries3 = fChain3->GetEntriesFast();
  Long64_t hgc_jentry = 0;

  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  Long64_t nbytes3 = 0, nb3 = 0;

  

  int decade = 0;
  int ahc_zeroHit = 0;



  bool DEBUG = false;
  //  bool DEBUG = true;
  int TOTAL_ACTIVE_LAYER = -1;
  int EE_LAYER = -1;
  int FH_LAYER = -1;
  int AH_LAYER = -1;
  if(!strcmp(conf_,"alpha") || !strcmp(conf_,"config1")) {
    TOTAL_ACTIVE_LAYER = 79;
    EE_LAYER = 28;
    FH_LAYER = 12;
    AH_LAYER = 39;
  }
  else if(!strcmp(conf_,"bravo") || !strcmp(conf_,"config2")){
    TOTAL_ACTIVE_LAYER = 78;
    EE_LAYER = 28;
    FH_LAYER = 11;
    AH_LAYER = 39;
  }
  else if(!strcmp(conf_,"charlie") || !strcmp(conf_,"config3")) {
    TOTAL_ACTIVE_LAYER = 59;
    EE_LAYER = 8;
    FH_LAYER = 12;
    AH_LAYER = 39;
    
  }
  else {
    cout<<"ERROR: Unknown configuration!!!!"<<endl;
    return;
  }

  //counter
  int nHgc=0, nAhc=0, nrechits=0;
  float offset = 159.4; //A constant to be added to AHC layer z positions
  //cout<<energy<<endl;
  int Min_ = atoi(folderID);
  //  int Max_ = atoi(max_);
  int chi2_method = atoi(chi_method);
  float FH_AH_relative_scale = 0.4;
  float alpha_ = FH_AH_relative_scale;
  float EE_scale = 94.624; //MIPs per Ge
  float FH_AH_scale = 12.788; //MIPs per GeV
  float ee_rescaling = 1.035;
  float fh_rescaling = 1.095;
  float ah_rescaling = 1.095;
  if(!strcmp(data,"data"))
    {
      ee_rescaling = 1;
      fh_rescaling = 1;
      ah_rescaling = 1;
    }

    cout<<"sim rescaling"<<"\t"<<ee_rescaling<<"\t"<<fh_rescaling<<"\t"<<ah_rescaling<<endl;
  char* outFileName = new char[1000];
  sprintf(outFileName,"./skimmed_ntuple_%s_chi2method%02d_%04d.root", data,chi2_method,Min_);//,Min_, Max_); // the name of the file where you write the tree. 
  TFile* outfile = TFile::Open(outFileName,"recreate");
  //TTree* pion_tree;
  pion_tree = new TTree("pion_variables_v1","variables for pion analysis v1");
  HGCNtupleVariables::init_piTree();
  int Elist[85] =  {10, 14 , 18 , 22 , 26 , 30,  34 , 38 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
                     82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
                     154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
                  226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
                     298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};
  char* name1 = new char[1000];
  //reading chi2 weights
   sprintf(name1,"./txt_maps/Enrange_min_max_SS_FH.txt");
   std::fstream file_H_;
   file_H_.open(name1,ios::in);
   sprintf(name1,"./txt_maps/Enrange_min_max_SS_EE.txt");
   std::fstream file_EH_;
   file_EH_.open(name1,ios::in);
   std::vector<map<int,double>> en_range_EE_xmin_vec;
   std::vector<map<int,double>> en_range_EE_xmax_vec;
   std::vector<map<int,double>> en_range_FH_xmin_vec;
   std::vector<map<int,double>> en_range_FH_xmax_vec;
   if(!file_H_.is_open())
        {
          std::cout << " file not opened" << std::endl;
        }
      else
        {
	  int energy;
          double min_, max_;
          while (file_H_ >> energy >>min_>>max_) {
            // int energy;
            // float min_, max_;
            // file_H_ >> energy >>min_>>max_;
            std::map<int,double> en_range_FH_xmin;
            std::map<int,double> en_range_FH_xmax;
            en_range_FH_xmin[energy]=min_;
            en_range_FH_xmax[energy]=max_;
            en_range_FH_xmin_vec.push_back(en_range_FH_xmin);
            en_range_FH_xmax_vec.push_back(en_range_FH_xmax);
          }
        }
   if(!file_EH_.is_open())
        {
          std::cout << " file not opened" << std::endl;
        }
      else
        {
	  int energy;
          double min_, max_;
          while (file_EH_ >>energy >>min_>>max_) {
            // int energy;
            // float min_, max_;
            // file_EH_ >>energy >>min_>>max_;
            std::map<int,double> en_range_FH_xmin;
            std::map<int,double> en_range_FH_xmax;
            en_range_FH_xmin[energy]=min_;
            en_range_FH_xmax[energy]=max_;
            en_range_EE_xmin_vec.push_back(en_range_FH_xmin);
            en_range_EE_xmax_vec.push_back(en_range_FH_xmax);
          }
        }

   sprintf(name1,"./parametres_weights_v1.txt");
   std::fstream file_p;
   file_p.open(name1,ios::in);
   std::vector<float> p0;
   std::vector<float> p1;
   std::vector<float> p2;
   if(!file_p.is_open())
     {
       std::cout << " file not opened" << std::endl;
     }
   else
     {
       float P0,P1,P2;
       while (file_p >>P0>>P1) {
	 p0.push_back(P0);
	 p1.push_back(P1);
	 //	 p2.push_back(P2);
       }
     }

   //fit function
   TF1* f_EH_trimAhcal_w1 = new TF1("f_EH_trimAhcal_w1","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
   f_EH_trimAhcal_w1->FixParameter(0,8.04233e-05);
   f_EH_trimAhcal_w1->FixParameter(1,1.62456);
   f_EH_trimAhcal_w1->FixParameter(2,1.31486);
   f_EH_trimAhcal_w1->FixParameter(3,1.14353);
   
  TF1* f_EH_trimAhcal_w2 = new TF1("f_EH_trimAhcal_w2","sqrt([0]*[0]+[1]*[1]/x)+exp(-[2]*x+[3])", 5, 355);
  f_EH_trimAhcal_w2->FixParameter(0,-3.49921e-07);
   f_EH_trimAhcal_w2->FixParameter(1,1.21772e+00);
   f_EH_trimAhcal_w2->FixParameter(2,1.72232e-04);
   f_EH_trimAhcal_w2->FixParameter(3,-3.50099e-02);
  
   TF1* f_EH_trimAhcal_w3= new TF1("f_EH_trimAhcal_w3","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
   f_EH_trimAhcal_w3->FixParameter(0,6.58726e-02);
   f_EH_trimAhcal_w3->FixParameter(1,2.62577);
   f_EH_trimAhcal_w3->FixParameter(2,-5.50631);
   f_EH_trimAhcal_w3->FixParameter(3,1.53569);

  TF1* f_H_trimAhcal_w1 = new TF1("f_H_trimAhcal_w1","sqrt([0]*[0]+[1]/x)", 5, 355);
  f_H_trimAhcal_w1->FixParameter(0,p0.at(6));
  f_H_trimAhcal_w1->FixParameter(1,p1.at(6));
  
  TF1* f_H_trimAhcal_w2= new TF1("f_H_trimAhcal_w2","sqrt([0]*[0]+[1]/x)+[2]/x", 5, 355);
  f_H_trimAhcal_w2->FixParameter(0,8.62610e-01);
  f_H_trimAhcal_w2->FixParameter(1,2.59516e+01);
  f_H_trimAhcal_w2->FixParameter(2,-7.02327e+00);
  TF1* f_H_trimAhcal_w3= new TF1("f_H_trimAhcal_w3","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
   f_H_trimAhcal_w3->FixParameter(0,2.05778e+00);
   f_H_trimAhcal_w3->FixParameter(1,1.12425e+01);
   f_H_trimAhcal_w3->FixParameter(2,-1.67933e+01);
   f_H_trimAhcal_w3->FixParameter(3,-4.91273e-01);


  // function for full ahcal scenarios
   TF1* f_EH_w1 = new TF1("f_EH_w1","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
   f_EH_w1->FixParameter(0,2.62110e-06);
   f_EH_w1->FixParameter(1,1.64811);
   f_EH_w1->FixParameter(2,1.18442);
   f_EH_w1->FixParameter(3,1.14035);

   TF1* f_EH_w2 = new TF1("f_EH_w2","sqrt([0]*[0]+[1]/x)+exp(-[2]*x+[3])", 5, 355);
   f_EH_w2->FixParameter(0,-3.15604e-06);//p0.at(4));
   f_EH_w2->FixParameter(1,1.27254e+00);//p1.at(4));
   f_EH_w2->FixParameter(2,1.55104e-04);//p2.at(4));
   f_EH_w2->FixParameter(3,-5.58071e-02);
   
   TF1* f_EH_w3 = new TF1("f_EH_w3","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
  f_EH_w3->FixParameter(0,-1.35815e-06);//p0.at(4));                                                                                       
  f_EH_w3->FixParameter(1,1.9127);//p1.at(4));                                                                                             
  f_EH_w3->FixParameter(2,-3.19441);//p2.at(4));                                                                                            
  f_EH_w3->FixParameter(3,1.01901);


   TF1* f_H_w1 = new TF1("f_H_w1","sqrt([0]*[0]+[1]/x)", 5, 355);
   f_H_w1->FixParameter(0,0);//p0.at());
   f_H_w1->FixParameter(1,0);//.at(6));
   TF1* f_H_w2= new TF1("f_H_w2","sqrt([0]*[0]+[1]/x)+[2]/x", 5, 355);
  f_H_w2->FixParameter(0,8.47109e-01);
  f_H_w2->FixParameter(1,2.38987e+01);
  f_H_w2->FixParameter(2,-6.49213e+00);
 

  TF1* f_H_w3= new TF1("f_H_w3","sqrt([0]*[0]+[1]*[1]/x)+[2]/x+[3]", 5, 355);
  f_H_w3->FixParameter(0,-6.04190e-07);
  f_H_w3->FixParameter(1,2.16114e+00);
  f_H_w3->FixParameter(2,-3.15370e+00);
  f_H_w3->FixParameter(3,9.17243e-01);

      TF1* fnct_EH_w1 = new TF1("fnct_EH_w1","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
      fnct_EH_w1->FixParameter(0,1.189897);
      fnct_EH_w1->FixParameter(1,5.215187);
      
      TF1* fnct_EH_w2 = new TF1("fnct_EH_w2","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
      fnct_EH_w2->FixParameter(0,1.004924);
      fnct_EH_w2->FixParameter(1,3.352913);
      TF1* fnct_EH_w3 = new TF1("fnct_EH_w3","sqrt([0]*[0]+[1]*[1]/x)", 5,355);
      fnct_EH_w3->FixParameter(0,1.025529);
      fnct_EH_w3->FixParameter(1,3.602263);
      
      TF1* fnct_H_w1 = new TF1("fnct_H_w1","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
      fnct_H_w1->FixParameter(0,0.0);
      fnct_H_w1->FixParameter(1,0.0);
      TF1* fnct_H_w2 = new TF1("fnct_H_w2","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
      fnct_H_w2->FixParameter(0,0.878231);
      fnct_H_w2->FixParameter(1,2.659706);
      TF1* fnct_H_w3 = new TF1("fnct_H_w3","sqrt([0]*[0]+[1]*[1]/x)", 5, 355);
      fnct_H_w3->FixParameter(0,0.984368);
      fnct_H_w3->FixParameter(1,2.765905);
    
  //  f_H_w2->FixParameter(2,p2.at(10));

  // TF1* f_H_w3 = new TF1("f_H_w3","sqrt([0]*[0]+[1]/x)", 5, 355);
  // f_H_w3->FixParameter(0, p0.at(11));
  // f_H_w3->FixParameter(1,p1.at(11));


  
   //declaring the matrix variables for chi2 optimization
  std::vector<ROOT::Math::SVector<double, 3> >consts_EH_vec;
  std::vector<ROOT::Math::SVector<double, 3> >values_EH_vec;
  std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > coeff_EH_vec;
  ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
  ROOT::Math::SVector<double, 3> consts_EH;
  ROOT::Math::SVector<double, 3> values_EH;

  bool isInverted_EH = false;
  for(int i =0; i<85; i++)
    {
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          coeffs_EH(i,j) = 0.0;
        }
        consts_EH(i) = 0.0;
        values_EH(i) = 0.0;
      }
        consts_EH_vec.push_back(consts_EH);
        values_EH_vec.push_back(values_EH);
        coeff_EH_vec.push_back(coeffs_EH);
    }
  std::vector<ROOT::Math::SVector<double, 2> >consts_H_vec;
  std::vector<ROOT::Math::SVector<double, 2> >values_H_vec;
  std::vector<ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2>> > coeff_H_vec;
  bool isInverted_H = false;
  ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
  ROOT::Math::SVector<double, 2> consts_H;
  ROOT::Math::SVector<double, 2> values_H;

  for(int i =0; i<85; i++)
    {
      for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
          coeffs_H(i,j) = 0.0;
        }
        consts_H(i) = 0.0;
        values_H(i) = 0.0;
      }
      consts_H_vec.push_back(consts_H);
      values_H_vec.push_back(values_H);
      coeff_H_vec.push_back(coeffs_H);
    }

  
  if(DEBUG) cout<<"DEBUG: Configuration = "<<conf_<<endl;
  if(DEBUG) cout<<"DEBUG: TOTAL_ACTIVE_LAYER = "<<TOTAL_ACTIVE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: EE_LAYER = "<<EE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: FH_LAYER = "<<FH_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: AH_LAYER = "<<AH_LAYER<<endl;

  if(DEBUG) cout<<"DEBUG: Entering event Loop"<<endl;

  Long64_t jentry;

  //reading the interaction length information layerwise
  float lambda[79];
  for(int i = 0; i < 79; i++) {
    lambda[i] = layer_positions[i+1].at(2); // for nuclear interaction length  & //pi_lambda use -at(3)
    //     cout<<lambda[i]<<endl;
  }

  // float lambda_pi[79];
  // for(int i = 0; i < 79; i++) {
    
  //   lambda_pi[i] = layer_positions[i+1].at(3); // for nuclear interaction length  & //pi_lambda use -at(3)
  //   cout<<lambda_pi[i]<<endl;
  // }
  int ahcal_layer[10]= {43,47,51,55,59,63,67,71,75,79}; 
  int En_bin[8]={20,50,80,100,120,200,250,300};
  for (jentry=0; jentry<nentries;jentry++,hgc_jentry++)
   {
         // for (jentry=0; jentry<1000;jentry++,hgc_jentry++) {
    // ==============print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;
    
      // ===============read this entry == == == == == == == == == == ==

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { break; cout<<"Breaking"<<endl;}
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout<<"insides the event loop: check1"<<endl;

      ////   MESSAGE ////
      // apparently code is not running beyond this point ///
      pi_event = 0;
      pi_run = 0;
      pi_pdgID = 0;
      pi_beamEnergy = -1;
      pi_trueBeamEnergy = -1;


      pi_energyLostEE = 0;
      pi_energyLostFH = 0;
      pi_energyLostBH = 0;
      pi_energyLostBeam = 0;
      pi_energyLostOutside=0;
      pi_energyLeakTransverseEE=0;
      pi_energyLeakTransverseFH=0;
      pi_energyLeakTransverseAH=0;
      pi_energyLeakLongitudinal=0;
      pi_energyLeakResidual=0;
      
      pi_recoMips=0;
      pi_recoMips_trimAhcal=0;
      pi_recoMips_AH_trimAhcal=0;
      pi_recoMips_EE=0;
      pi_recoMips_FH=0;
      pi_recoMips_AH=0;
      
      pi_recoFixwt=0;
      pi_recoFixwt_trimAhcal=0;
      pi_recoFixwt_EE=0;
      pi_recoFixwt_FH=0;
      pi_recoFixwt_AH=0;
      pi_recoFixwt_AH_trimAhcal=0;
      pi_chi2Reco=0;
      pi_chi2Reco_EE=0;
      pi_chi2Reco_FH=0;
      pi_chi2Reco_AH=0;
      pi_trimAhcal_chi2Reco=0;
      pi_trimAhcal_chi2Reco_EE=0;
      pi_trimAhcal_chi2Reco_FH=0;
      pi_trimAhcal_chi2Reco_AH=0;
      pi_rechit_shower_start_layer=-1;
      
      event_count[1]++;
      // if(!pi_isHGC_AHC_sync)
      // 	{ cout<<"wrongevents"<<endl;}
      //if(Nrechits_FH_module_42 >80 ||  Nrechits_FH_module_45 > 80){count_fh++; continue;}
      if(isFHNoisy) {count_fh++; continue;}
      event_count[2]++;
      h_true_beamenergy[1]->Fill(trueBeamEnergy);
       //good track      
      bool isGoodTrack = (dwcReferenceType >= 13 && trackChi2_X < 10 && trackChi2_Y < 10);
      if(!isGoodTrack) { count_badtrack++; continue;}
      event_count[3]++;
      h_true_beamenergy[2]->Fill(trueBeamEnergy);
      //      muon veto
      if(MuonVeto) continue;
      event_count[4]++;
      h_true_beamenergy[3]->Fill(trueBeamEnergy);
      //preshower
      if(rechit_shower_start_layer<=2) continue;
      event_count[5]++;
      h_true_beamenergy[4]->Fill(trueBeamEnergy);
      // //track impact only for TB data & discrete simulation but not for flat energy samples (no track impact information was available)
      if(!strcmp(data,"data") || !strcmp(data,"sim_discrete")) {
      	if(!isInTrackWindow) continue;
      }
      event_count[6]++;
      h_beamenergy->Fill(beamEnergy);
      h_particle->Fill(pdgID);
      
      //for tree
      //pi_isHGC_AHC_sync = true;
      pi_event = event;
      pi_run = run;
      pi_pdgID = pdgID;
      
      pi_beamEnergy = beamEnergy;
      pi_trueBeamEnergy = trueBeamEnergy;
      //pi_NRechits = NRechits;
      // pi_rechit_cogX=rechit_cogX;
      // pi_rechit_cogY=rechit_cogY;
      // pi_impactX_layer=TrackImpactX_layer;
      // pi_impactY_layer=TrackImpactY_layer;
      // pi_rechit_energyPerLayer=rechit_energyPerLayer;
      // pi_rechit_nHitsPerLayer=rechit_nHitsPerLayer;
      pi_rechit_shower_start_layer=rechit_shower_start_layer;
      // pi_ntracks = ntracks;
      // pi_trackChi2_X = trackChi2_X;
      // pi_trackChi2_Y = trackChi2_Y;
      // pi_dwcReferenceType = dwcReferenceType;
      // pi_m_x = m_x;
      // pi_m_y = m_y;
      // pi_b_x = b_x;
      // pi_b_y = b_y;

      //pi_ahc_nHits = ahc_nHits;
      // pi_ahc_energyPerLayer = ahc_energyPerLayer;
      // pi_ahc_nHitsPerLayer = ahc_nHitsPerLayer;
      pi_energyLostEE = energyLostEE;
      pi_energyLostFH = energyLostFH;
      pi_energyLostBH = energyLostBH;
      pi_energyLostBeam = energyLostBeam;
      
      //pion_tree->Fill();

      //variables declaration
      totalEnergy_inGeV=0;
      EnergySum_SSinEE=0;
      EnergySum_SSinFH=0;
      Esum_rechits_FH=0;
      Esum_rechits_EE=0;
      Esum_rechits_AH=0;
      Esum_rechits_FH_inGeV=0;
      Esum_rechits_EE_inGeV=0;
      Esum_rechits_AH_inGeV=0;
      
      ////////////////////////////////////////////
      //            HGCAL Part                  //
      ////////////////////////////////////////////
      total_rechits=0;
      rechits_EE=0;
      rechits_FH=0;
      rechits_AH=0;
      int nrechit_trimAhcal=0;
      /// Read HGCAL Tree
      if(DEBUG) cout<<"DEBUG: Start Analylizing HGCAL  RecHits!!"<<endl;
      // total_rechits_energy_EE=0;
      // total_rechits_energy_FH=0;
      float dr_a=999.0,dr_b=999.0;
      int count_rechit[4]={};
      for( int i = 0 ; i < NRechits; i++)
      	{
      	  int channel = rechit_channel->at(i);
      	  int chip = rechit_chip->at(i);
      	  int en_chan = 1000*chip + channel;
      	  int layer = rechit_layer->at(i);
	  float energy=rechit_energy->at(i);
	  if(hgc_channel_mask->at(i)) continue;
  	  if(!pass_noise_thres->at(i))continue;
	  double trackx = TrackImpactX_layer->at(layer-1);
          double tracky =TrackImpactY_layer->at(layer-1);
	  float lam_ = lambda[layer-1];
	  count_rechit[0]++;
	  // pi_lambda.push_back(lam_);
	  // pi_lambda_trimAhcal.push_back(lam_);
	  float recx = rechit_x->at(i);
	  float	recy = rechit_y->at(i);
	  float recx_fin = recx;
          float recy_fin = recy;
	  // regularizing the hits going beyond HGCAl physical boundaries
	  
	  // pi_rechit_layer.push_back(layer);
	  // pi_rechit_energy.push_back(energy);
  	  // pi_rechit_x.push_back(recx);//hit_x->at(i));
  	  // pi_rechit_y.push_back(recy);//hit_y->at(i));
	  // pi_rechitEn_trimAhcal.push_back(energy);
	  // pi_comb_rechit_x_trimAhcal.push_back(recx);
	  // pi_comb_rechit_y_trimAhcal.push_back(recy);
	  // pi_comb_rechit_z_trimAhcal.push_back(rechit_z->at(i));
  	  // pi_rechit_z.push_back(rechit_z->at(i));
  	  // pi_rechit_iu.push_back(rechit_iu->at(i));
  	  // pi_rechit_iv.push_back(rechit_iv->at(i));
  	  // pi_rechit_iU.push_back(rechit_iU->at(i));
  	  // pi_rechit_iV.push_back(rechit_iV->at(i));
  	  // pi_rechit_amplitudeHigh.push_back(rechit_amplitudeHigh->at(i));
  	  // pi_rechit_amplitudeLow.push_back(rechit_amplitudeLow->at(i));
  	  // pi_rechit_noise_flag.push_back(rechit_noise_flag->at(i));
	  // pi_rechit_module.push_back(rechit_module->at(i));
  	  // pi_rechit_chip.push_back(rechit_chip->at(i));
  	  // pi_rechit_channel.push_back(rechit_channel->at(i));
  	  // pi_rechit_type.push_back(rechit_type->at(i));

	  // pi_combined_rechits_energy.push_back(energy);
	  // pi_combined_rechit_x.push_back(recx);//hit_x->at(i));
	  // pi_combined_rechit_y.push_back(recy);//hit_y->at(i));
	  // pi_combined_rechit_z.push_back(rechit_z->at(i));
	  nHgc++;
	  nrechits++;
	  nrechit_trimAhcal++;
	  
	  if(layer<=28)
	    {
	      
	      energy = energy/ee_rescaling;
	      h_rechits_mips_EE->Fill(energy);
	      h_rechits_GeV_EE->Fill(energy*0.0105);
	      Esum_rechits_EE+=energy;
	      rechits_EE++;
	    }
	  else
	    {
	      energy = energy/fh_rescaling;
	      h_rechits_mips_FH->Fill(energy);
              h_rechits_GeV_FH->Fill(energy*0.0789);
	      Esum_rechits_FH+=energy;
	      rechits_FH++;
	    }
  	  } //Nrechits loop
      // pi_NRechits = nHgc;

      
      // ////////////////////////////////////////////
      // //            AHCAL Part                  //
      // ////////////////////////////////////////////
      
      /// Read AHCAL Tree
      
      if(DEBUG) cout<<"DEBUG: Start Analylizing AHCAL  RecHits!!"<<endl;
      //cout<<ahc_nHits<<"ahcal"<<endl;
      float rechitEnergySum_AH_trimm=0.0;
      for(int i = 0 ; i < ahc_nHits; i++) {
      	if(ahc_channel_mask->at(i)) continue;
  	// pi_ahc_hitI.push_back(ahc_hitI->at(i));
  	// pi_ahc_hitJ.push_back(ahc_hitJ->at(i));
  	// pi_ahc_hitK.push_back(ahc_hitK->at(i));
	float lam_ = lambda[ahc_hitK->at(i)+39];
	// float lam_pi = lambda_pi[ahc_hitK->at(i)+39];
	// pi_lambda_pi.push_back(lam_pi);
	// pi_lambda.push_back(lam_);
	count_rechit[0]++;
 	// pi_ahc_hitZ.push_back(ahc_hitZ->at(i));
  	// pi_ahc_hitEnergy.push_back(ahc_hitEnergy->at(i));
	int layer = ahc_hitK->at(i)+40;
	//alignment corrections
	double trackx = TrackImpactX_layer->at(layer-1);// track_x[layer-1];
	double tracky =  TrackImpactY_layer->at(layer-1);//track_y[layer-1];
	
	float recx =ahc_hitX->at(i);
	float recy =ahc_hitY->at(i);
	// pi_ahc_hitX.push_back(recx);
        // pi_ahc_hitY.push_back(recy);

	// pi_combined_rechits_energy.push_back(ahc_hitEnergy->at(i));
	// pi_combined_rechit_x.push_back(recx);
	// pi_combined_rechit_y.push_back(recy);
	// pi_combined_rechit_z.push_back(ahc_hitZ->at(i)+offset);
	nAhc++;
	nrechits++;

	Esum_rechits_AH+=(ahc_hitEnergy->at(i)/ah_rescaling);
	h_rechits_mips_AH->Fill(ahc_hitEnergy->at(i));
	
	h_rechits_GeV_AH->Fill(ahc_hitEnergy->at(i)*0.0316);

	int count=0;
	for(int i_l=0; i_l<10; i_l++)
	  {
	    // if(count==1)
	    //   break;
	    if(layer ==ahcal_layer[i_l])
	      {
		// pi_trimmed_ahcal_flag.push_back(1);
		//pi_lambda_pi_trimAhcal.push_back(lam_pi);
		// pi_lambda_trimAhcal.push_back(lam_);
		// pi_rechitEn_trimAhcal.push_back(ahc_hitEnergy->at(i));
		// pi_comb_rechit_x_trimAhcal.push_back(recx);
		// pi_comb_rechit_y_trimAhcal.push_back(recy);
		// pi_comb_rechit_z_trimAhcal.push_back(ahc_hitZ->at(i)+offset);
		rechitEnergySum_AH_trimm+=(ahc_hitEnergy->at(i)/ah_rescaling);
		nrechit_trimAhcal++;
		count++;		
	      }
	  }
	//cout<<count<<endl;
	//if(count==0)
	//	pi_trimmed_ahcal_flag.push_back(count);
	  
      }// AHCAL hit loop
   
      // code to get the updated relative weight between FH and AH after emoving additional layers
      h_true_beamenergy[5]->Fill(trueBeamEnergy);
      // pi_ahc_nHits=nAhc;
      // pi_combined_NRechits=nrechits;
      // pi_Nrechit_trimAhcal=nrechit_trimAhcal;
      nrechit_trimAhcal=0;
      //      pion_tree->Fill();
      // cout<<"total rechits"<<"\t"<<count_rechit[0]<<"\t"<<count_rechit[1]<<"\t"<<count_rechit[2]<<endl;
      // cout<<"fraction"<<"\t"<<count_rechit[1]/count_rechit[0]<<"\t"<<count_rechit[2]/count_rechit[0]<<endl;
      
      count_rechit[0]=0;
      count_rechit[1]=0;
      count_rechit[2]=0;
      count_rechit[3]=0;
      rechits_AH = nAhc;
      total_rechits = nrechits;
      nrechits=0;
      nHgc=0;
      nAhc =0;
      //enegry in MIPs
      pi_recoMips=Esum_rechits_EE+Esum_rechits_FH+Esum_rechits_AH;
      pi_recoMips_trimAhcal=Esum_rechits_EE+Esum_rechits_FH+rechitEnergySum_AH_trimm;
      pi_recoMips_AH_trimAhcal=rechitEnergySum_AH_trimm;
      pi_recoMips_EE=Esum_rechits_EE;
      pi_recoMips_FH=Esum_rechits_FH;
      pi_recoMips_AH=Esum_rechits_AH;
      
      Esum_rechits_EE_inGeV=0.0105*Esum_rechits_EE;
      Esum_rechits_FH_inGeV= 0.0789*Esum_rechits_FH;
      Esum_rechits_AH_inGeV=0.0316*Esum_rechits_AH;
      //rechit in gev with fixwt
      pi_recoFixwt = Esum_rechits_EE_inGeV + Esum_rechits_FH_inGeV + Esum_rechits_AH_inGeV;
      pi_recoFixwt_EE = Esum_rechits_EE_inGeV;
      pi_recoFixwt_FH =	Esum_rechits_FH_inGeV;
      pi_recoFixwt_AH =	Esum_rechits_AH_inGeV;
      float Esum_rechits_trimAH_inGeV = 0.0757*rechitEnergySum_AH_trimm;
      pi_recoFixwt_trimAhcal = Esum_rechits_EE_inGeV + Esum_rechits_FH_inGeV + Esum_rechits_trimAH_inGeV;
      pi_recoFixwt_AH_trimAhcal = Esum_rechits_trimAH_inGeV;
      
      double rechitEnergySum_EE = Esum_rechits_EE;
      double rechitEnergySum_FH = Esum_rechits_FH;
      double rechitEnergySum_AH = Esum_rechits_AH;
      double EE_detscale = Esum_rechits_EE_inGeV;//(0.0105rechitEnergySum_EE/EE_scale);                                        
      double FH_detscale = Esum_rechits_FH_inGeV;//(rechitEnergySum_FH/FH_AH_scale);                                                 
      double AH_detscale = Esum_rechits_AH_inGeV;//(alpha_*rechitEnergySum_AH)/FH_AH_scale;
      double trimAhcal_detscale = Esum_rechits_trimAH_inGeV;
      double full_energy = FH_detscale+AH_detscale;                                                                                       
      double total_energy= EE_detscale + FH_detscale + AH_detscale;//Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_trimAH_inGeV;//(rechitEnergySum_EE/EE_scale) + (rechitEnergySum_FH+ (alpha_ * rechitEnergySum_AH))/FH_AH_scale;
      double tot_E_gev= total_energy;
      double total_energy_trim = EE_detscale + FH_detscale+trimAhcal_detscale;
      ///////////////////////////////////////////////////////////                                                                           
      /////     H hadrons ; Chi2 matrix initialzation     /////                                                                             
      //////////////////////////////////////////////////////////
      //filling the chi2 energy
      //cout<<"check2"<<"H-hadronsinitialize tree branch"<<endl;
      //total_energy = (rechitEnergySum_FH+ (alpha_ * rechitEnergySum_AH))/FH_AH_scale+0.4;
      // total_energy =total_energy_trim;
      // AH_detscale = trimAhcal_detscale;
      // tot_E_gev= total_energy;
      double E_beam =0.0;
      double O = 0.0;
      double sigma2_H = (0.089*0.089)*(total_energy*total_energy) + (1.228*1.228)*total_energy;
      //double tot_E_gev= total_energy;
      
      //      cout<<rechit_shower_start_layer<<"\t"<<trueBeamEnergy<<endl;
      for(int i_bin=0;i_bin<85;i_bin++)
	{
	  //cout<<"check4"<<"in the loop"<<"\t"<<i_bin<<endl;
	  bool IsIn_en_range = false;
          if(trueBeamEnergy>=Elist[i_bin] && trueBeamEnergy <=Elist[i_bin]+4)
            {
	      //cout<<"check7"<<"in the loop"<<"\t"<<i_bin<<endl;
	      
	      E_beam = Elist[i_bin]+2;
	      
	      if(rechit_shower_start_layer>28)
		{

		  float w1=0,w2=0,w3=0;
		  //cout<<"check8"<<"in the loop"<<"\t"<<i_bin<<endl;
		  if(!strcmp(data,"sim_discrete"))
		    {
		       w1 = fnct_H_w1->Eval(trueBeamEnergy);
		       w2 = fnct_H_w2->Eval(trueBeamEnergy);
		       w3 = fnct_H_w3->Eval(trueBeamEnergy);

		    }
		  else if (!strcmp(data,"data"))
		    {
		      trueBeamEnergy = beamEnergy;
		    
		       w1 = f_H_w1->Eval(trueBeamEnergy);
		       w2 = f_H_w2->Eval(trueBeamEnergy);
		       w3 = f_H_w3->Eval(trueBeamEnergy);
		    }
		  else
		    {
		       w1 = f_H_w1->Eval(trueBeamEnergy);
                       w2 = f_H_w2->Eval(trueBeamEnergy);
                       w3 = f_H_w3->Eval(trueBeamEnergy);

		    }
		  float chi2_funct = EE_detscale +w2*FH_detscale + w3*AH_detscale;
		  hist_resp_total_funct[i_bin]->Fill(chi2_funct);
		  hist_resp_SS_FH_funct[i_bin]->Fill(chi2_funct);
		  pi_chi2Reco = chi2_funct;
		  pi_chi2Reco_EE = EE_detscale;
		  pi_chi2Reco_FH = w2*FH_detscale;
		  pi_chi2Reco_AH = w3*AH_detscale;
		  if (!strcmp(data,"data"))
                    
                      trueBeamEnergy = beamEnergy;

		  float trmiAhcal_w1 = f_H_trimAhcal_w1->Eval(trueBeamEnergy);
		  float trimAhcal_w2 = f_H_trimAhcal_w2->Eval(trueBeamEnergy);
		  float trimAhcal_w3 = f_H_trimAhcal_w3->Eval(trueBeamEnergy);
		  chi2_funct = EE_detscale +trimAhcal_w2*FH_detscale + trimAhcal_w3*trimAhcal_detscale;
		  hist_resp_total_trimAhcal[i_bin]->Fill(chi2_funct);
		  hist_resp_SS_FH_trimAhcal[i_bin]->Fill(chi2_funct);
		  pi_trimAhcal_chi2Reco = chi2_funct;
                  pi_trimAhcal_chi2Reco_EE = EE_detscale;
                  pi_trimAhcal_chi2Reco_FH = trimAhcal_w2*FH_detscale;
                  pi_trimAhcal_chi2Reco_AH = trimAhcal_w3*trimAhcal_detscale;

		   w1 = getChi2Weights_H(E_beam).at(0);
		   w2 = getChi2Weights_H(E_beam).at(1);
		   w3 = getChi2Weights_H(E_beam).at(2);
		  
		  //cout<<"check5"<<"in the loop"<<"\t"<<i_bin<<endl;
		  float chi2_energy = rechitEnergySum_EE + w2*rechitEnergySum_FH + w3*rechitEnergySum_AH ;
		  if(chi2_method > 0)  {
		    if(chi2_method == 1) O = 0.4;
		    chi2_energy = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale ;
		  }
		  hist_resp_total_1[i_bin]->Fill(chi2_energy);
		  hist_resp_SS_FH_1[i_bin]->Fill(chi2_energy);
		  hist_fixwt_total[i_bin]->Fill(total_energy);
                  hist_fixwt_SS_FH[i_bin]->Fill(total_energy);

		  // ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
		  // ROOT::Math::SVector<double, 2> consts_H;
		  // ROOT::Math::SVector<double, 2> values_H;
		  // for(int i = 0; i < 2; i++) {
		  //   for(int j = 0; j < 2; j++) {
		  //     coeffs_H(i,j) = 0.0;
		  //   }	
		  //   consts_H(i) = 0.0;
		  //   values_H(i) = 0.0;
		  // }
		  // IsIn_en_range = (total_energy > en_range_FH_xmin_vec[i_bin][Elist[i_bin]+2] && total_energy < en_range_FH_xmax_vec[i_bin][Elist[i_bin]+2]);
		  // //IsIn_en_range = true;
		  // if(IsIn_en_range){

		  // if(chi2_method == 0)
		  //   {
		  //     coeffs_H(0,0) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2_H;
		  //     coeffs_H(0,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2_H;
		  //     coeffs_H(1,0) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2_H;
		  //     coeffs_H(1,1) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2_H;
		      
		  //     consts_H(0)   = (E_beam-O)*rechitEnergySum_FH/sigma2_H;
		  //     consts_H(1)   = (E_beam-O)*rechitEnergySum_AH/sigma2_H;
		  //     consts_H_vec[i_bin]+= consts_H;
		  //     coeff_H_vec[i_bin]+=coeffs_H;
		  //   }
		  // else
		  //   {
		  //     coeffs_H(0,0) = (FH_detscale*FH_detscale)/sigma2_H;
		  //     coeffs_H(0,1) = (AH_detscale*FH_detscale)/sigma2_H;
		  //     coeffs_H(1,0) = (FH_detscale*AH_detscale)/sigma2_H;
		  //     coeffs_H(1,1) = (AH_detscale*AH_detscale)/sigma2_H;
		      
		  //     consts_H(0)   = (E_beam-O)*FH_detscale/sigma2_H;
		  //     consts_H(1)   = (E_beam-O)*AH_detscale/sigma2_H;
		  //     consts_H_vec[i_bin] += consts_H;
		  //     coeff_H_vec[i_bin]+=coeffs_H;

		  //   }
		}
		
	      else if (rechit_shower_start_layer<=28)
		{
		  float w1_func=0.0, w2_func=0.0, w3_func=0.0;
		  //cout<<"check6"<<"in the loop"<<"\t"<<i_bin<<endl;
		  if(!strcmp(data,"sim_discrete"))
                    {
                       w1_func = fnct_H_w1->Eval(trueBeamEnergy);
                       w2_func = fnct_H_w2->Eval(trueBeamEnergy);
                       w3_func = fnct_H_w3->Eval(trueBeamEnergy);

                    }
                  else if (!strcmp(data,"data"))
                    {
                      trueBeamEnergy = beamEnergy;

                       w1_func = f_H_w1->Eval(trueBeamEnergy);
                       w2_func = f_H_w2->Eval(trueBeamEnergy);
                       w3_func = f_H_w3->Eval(trueBeamEnergy);
                    }
		  else
		    {
		       w1_func = f_EH_w1->Eval(trueBeamEnergy);
		       w2_func = f_EH_w2->Eval(trueBeamEnergy);
		       w3_func = f_EH_w3->Eval(trueBeamEnergy);
		    }
		  float chi2_funct = w1_func*EE_detscale +w2_func*FH_detscale + w3_func*AH_detscale;
                  hist_resp_total_funct[i_bin]->Fill(chi2_funct);
                  hist_resp_SS_EE_funct[i_bin]->Fill(chi2_funct);
		  pi_chi2Reco =	chi2_funct;
		  pi_chi2Reco_EE = w1_func*EE_detscale;
                  pi_chi2Reco_FH = w2_func*FH_detscale;
                  pi_chi2Reco_AH = w3_func*AH_detscale;
		  if (!strcmp(data,"data"))
                  
		    trueBeamEnergy = beamEnergy;

                  
                  float trimAhcal_w1 = f_EH_trimAhcal_w1->Eval(trueBeamEnergy);
                  float trimAhcal_w2 = f_EH_trimAhcal_w2->Eval(trueBeamEnergy);
                  float trimAhcal_w3 = f_EH_trimAhcal_w3->Eval(trueBeamEnergy);
		  chi2_funct = trimAhcal_w1*EE_detscale +trimAhcal_w2*FH_detscale + trimAhcal_w3*trimAhcal_detscale;
                  hist_resp_total_trimAhcal[i_bin]->Fill(chi2_funct);
                  hist_resp_SS_EE_trimAhcal[i_bin]->Fill(chi2_funct);
		  pi_trimAhcal_chi2Reco = chi2_funct;
                  pi_trimAhcal_chi2Reco_EE = trimAhcal_w1*EE_detscale;
                  pi_trimAhcal_chi2Reco_FH = trimAhcal_w2*FH_detscale;
                  pi_trimAhcal_chi2Reco_AH = trimAhcal_w3*trimAhcal_detscale;

		  float w1 = getChi2Weights_EH(E_beam).at(0);                                                                          
		  float w2 = getChi2Weights_EH(E_beam).at(1);                                                                         
		  float w3 = getChi2Weights_EH(E_beam).at(2);
		  //cout<<"check6"<<"in the loop"<<"\t"<<i_bin<<endl;
		  float chi2_energy = w1*EE_detscale + w2*EE_detscale + w3*EE_detscale;// + O;
		  if(chi2_method > 0)  {
                    if(chi2_method == 1) O = 0.4;
                    chi2_energy = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale;// + O;
                  }
		  
		  hist_resp_total_1[i_bin]->Fill(chi2_energy);
                  hist_resp_SS_EE_1[i_bin]->Fill(chi2_energy);
		  hist_fixwt_total[i_bin]->Fill(tot_E_gev);
                  hist_fixwt_SS_EE[i_bin]->Fill(tot_E_gev);


		  // double sigma2_EH =(0.084*0.084)*(tot_E_gev*tot_E_gev) + (1.39*1.39)*tot_E_gev;// (0.084*0.084) + (1.39*1.39)*tot_E_gev;
		  // ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
		  // ROOT::Math::SVector<double, 3> consts_EH;
		  // ROOT::Math::SVector<double, 3> values_EH;

		  // for(int i = 0; i < 3; i++) {
		  //   for(int j = 0; j < 3; j++) {
		  //     coeffs_EH(i,j) = 0.0;
		  //   }
		  //   consts_EH(i) = 0.0;
		  //   values_EH(i) = 0.0;
		  // }
		  
		  //  IsIn_en_range = (tot_E_gev > en_range_EE_xmin_vec[i_bin][Elist[i_bin]+2] && tot_E_gev < en_range_EE_xmax_vec[i_bin][Elist[i_bin]+2]);
		  //  //IsIn_en_range = true;
		  //  if(IsIn_en_range){
		  //    if(chi2_method==0)
		  //      {
		  //     coeffs_EH(0,0) = (rechitEnergySum_EE*rechitEnergySum_EE)/sigma2_EH;
		  //     coeffs_EH(0,1) = (rechitEnergySum_EE*rechitEnergySum_FH)/sigma2_EH;
		  //     coeffs_EH(0,2) = (rechitEnergySum_EE*rechitEnergySum_AH)/sigma2_EH;
		  //     coeffs_EH(1,0) = (rechitEnergySum_FH*rechitEnergySum_EE)/sigma2_EH;
		  //     coeffs_EH(1,1) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2_EH;
		  //     coeffs_EH(1,2) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2_EH;
		  //     coeffs_EH(2,0) = (rechitEnergySum_AH*rechitEnergySum_EE)/sigma2_EH;
		  //     coeffs_EH(2,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2_EH;
		  //     coeffs_EH(2,2) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2_EH;

		  //     consts_EH(0)   = E_beam*rechitEnergySum_EE/sigma2_EH;
		  //     consts_EH(1)   = E_beam*rechitEnergySum_FH/sigma2_EH;
		  //     consts_EH(2)   = E_beam*rechitEnergySum_AH/sigma2_EH;
		  //     consts_EH_vec[i_bin]+= consts_EH;
		  //     coeff_EH_vec[i_bin]+=coeffs_EH;

		  //   }
		  // else
		  //   {
		  //     coeffs_EH(0,0) = (EE_detscale*EE_detscale)/sigma2_EH;
		  //     coeffs_EH(0,1) = (EE_detscale*FH_detscale)/sigma2_EH;
		  //     coeffs_EH(0,2) = (EE_detscale*AH_detscale)/sigma2_EH;
		  //     coeffs_EH(1,0) = (FH_detscale*EE_detscale)/sigma2_EH;
		  //     coeffs_EH(1,1) = (FH_detscale*FH_detscale)/sigma2_EH;
		  //     coeffs_EH(1,2) = (FH_detscale*AH_detscale)/sigma2_EH;
		  //     coeffs_EH(2,0) = (AH_detscale*EE_detscale)/sigma2_EH;
		  //     coeffs_EH(2,1) = (AH_detscale*FH_detscale)/sigma2_EH;
		  //     coeffs_EH(2,2) = (AH_detscale*AH_detscale)/sigma2_EH;


		  //     consts_EH(0)   = (E_beam - O)*EE_detscale/sigma2_EH;
		  //     consts_EH(1)   = (E_beam - O)*FH_detscale/sigma2_EH;
		  //     consts_EH(2)   = (E_beam - O)*AH_detscale/sigma2_EH;
		  //     consts_EH_vec[i_bin]+= consts_EH;
		  //     coeff_EH_vec[i_bin]+=coeffs_EH;
		      
		  //   }
		  //  }
		}
            }	         
        }//true bin loop
      //cout<<"check3"<<"After chi2"<<endl;
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      
      //if(DEBUG && jentry > 10) break;
      //      gSystem->Exit(0);
      //if(jentry > 2) break;
      pion_tree->Fill();
    } // loop over entries
  
//   char* name = new char[1000];
//   sprintf(name,"chi2_calibFact_H_hadrons_UpdatedflatEn_%d_scalMC_2sigma_offSetEE_DownscaledAHCAL_v1.txt",chi2_method);
//   std::ofstream file_H;
//   file_H.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_H_hadrons_flatEn_%d.txt",chi2_method,ios::out);                        
//   sprintf(name,"chi2_calibFact_EH_hadrons_UpdatedflatEn_%d_scalMC_2sigma_offSetEE_DownscaledAHCAL_v1.txt",chi2_method);
//   std::ofstream file_EH;
//   file_EH.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_EH_hadrons_flatEn_%d.txt",chi2_method,ios::out

//   /////////////////////////////////////////////////////                                                                                     
//   //////                 EH Hadrons            ////////                                                                                     
//   /////////////////////////////////////////////////////                                                                                     

//    for(int i_en =0; i_en<85; i_en++)
//     {
//       ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
//       ROOT::Math::SVector<double, 3> consts_EH;
//       ROOT::Math::SVector<double, 3> values_EH;
//       ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
//       ROOT::Math::SVector<double, 2> consts_H;
//       ROOT::Math::SVector<double, 2> values_H;
//        //      bool isInverted_EH = false;                                                                                                  
//       for(int i = 0; i < 3; i++) {
//         for(int j = 0; j < 3; j++) {
//           coeffs_EH(i,j) = 0.0;
//         }
//         consts_EH(i) = 0.0;
//         values_EH(i) = 0.0;
//       }
//        for(int i = 0; i < 2; i++) {
//         for(int j = 0; j < 2; j++) {
//           coeffs_H(i,j) = 0.0;
//         }
//         consts_H(i) = 0.0;
//         values_H(i) = 0.0;
//       }

//       consts_EH = consts_EH_vec[i_en];
//       values_EH = values_EH_vec[i_en];
//       coeffs_EH = coeff_EH_vec[i_en];
//       isInverted_EH = coeffs_EH.Invert();
//       for(int i = 0; i < 3; i++) {
//         for(int j = 0; j < 3; j++) {
//           cout<<coeffs_EH(i, j)<<"\t";
//         }
//         cout<<endl;
//       }

//       cout<<" "<<consts_EH(0)<<"\t"<<consts_EH(1)<<"\t"<<consts_EH(2)<<endl;
// consts_H = consts_H_vec[i_en];
//       values_H = values_H_vec[i_en];
//       coeffs_H = coeff_H_vec[i_en];

//       isInverted_H = coeffs_H.Invert();
//       for(int i = 0; i <3; i++) {
//         for(int j = 0; j < 3; j++) {
//           cout<<coeffs_H(i, j)<<"\t";
//         }
//         cout<<endl;
//       }
//       cout<<" "<<consts_H(0)<<"\t"<<consts_H(1)<<"\t"<<consts_H(2)<<endl;
//       if(isInverted_EH) {
//         values_EH = coeffs_EH*consts_EH;
//         cout<<"EH Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = "<<values_EH(0)<<" ;w2 = "<<values_EH(1)<<" ;w3 = "<<values_EH(2)<<";"<<endl;
// 	cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_EH(0,0))<< " " << sqrt(coeffs_EH(1,1)) << " " << sqrt(coeffs_EH(2,2))<<endl;

//          file_EH<<Elist[i_en]+2<<"\t"<<values_EH(0)<<"\t"<<values_EH(1)<<"\t"<<values_EH(2)<<"\t"<<sqrt(coeffs_EH(0,0))<<"\t"<<sqrt(coeffs_EH(1,1))<<"\t"<<sqrt(coeffs_EH(2,2))<<"\n";

//         //file_EH<<Elist[i_en]+2<<"\t"<<values_EH(0)<<"\t"<<values_EH(1)<<"\t"<<values_EH(2)<<"\n";
//       }
//       else {
//         cout<<"Error: Could not invert for EH hadrons..."<<endl;
//       }
//       //isInverted_H = true;                                                                                                                
//       if(isInverted_H) {
//         values_H = coeffs_H*consts_H;
//         //cout<<"H Hadrons => For E = "<<E_beam<<"GeV, w1 = "<<values_H(0)<<" ;w2 = "<<values_H(1)<<" ;w3 = "<<values_H(2)<<";"<<endl;      
//         cout<<"H Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = 0.0 ;w2 = "<<values_H(0)<<" ;w3 = "<<values_H(1)<<";"<<endl;
// 	cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_H(0,0)) << " " <<sqrt(coeffs_H(1,1))<<endl;

//         file_H<<Elist[i_en]+2<<"\t"<<0<<"\t"<<values_H(0)<<"\t"<<values_H(1)<<"\t"<<0<<"\t"<<sqrt(coeffs_H(0,0))<<"\t"<<sqrt(coeffs_H(1,1))<<"\n";

//         //file_H<<Elist[i_en]+2<<"\t"<<values_H(0)<<"\t"<<values_H(1)<<"\n";
//       }
//       else {
//         cout<<"Error: Could not invert for H hadrons..."<<endl;
//       }
//     }


  cout<<endl<<endl;


  ///////////////////////////////////////////////////////////
  ///////  E N D     O F     E N T R Y     L O O P     //////
  ///////////////////////////////////////////////////////////
  //  gSystem->Exit(0);
  cout<<"Got Out "<<jentry<<endl;
  cout<<count_fh<<endl;
  cout<<count_badtrack<<endl;
  for(int i =0;i<7;i++)
    {
      cout<<event_count[i]<<endl;
     }
  outfile->cd();
  pion_tree->Write();
  outfile->Close();
  cout<<"outFile: "<<outFileName<<" written!!"<<endl;
}
