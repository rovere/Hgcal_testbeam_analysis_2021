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

  if (argc < 7) {
    cerr << "Please give 7 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" " << "configuration" 
	 <<" " << "energy" << "min_"<<"max_"<< endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *config          = argv[4];
  const char *energy = argv[5];
  const char *min_ = argv[6];
  const char *max_ = argv[7];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy,min_,max_);
  cout << "dataset " << data << " " << endl;
  cout << "configuration " << config << " " << endl;
  cout << "energy " << energy << " " << endl;
 int chi2_method = atoi(energy);
  hgcOctTB.EventLoop(data,energy, min_, max_);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data, const char *energy, const char *min_, const char *max_) {


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
  int Min_ = atoi(min_);
  int Max_ = atoi(max_);
  int chi2_method = atoi(energy);
  float FH_AH_relative_scale = 0.4;
  float alpha_ = FH_AH_relative_scale;
  float EE_scale = 94.624; //MIPs per Ge
  float FH_AH_scale = 12.788; //MIPs per GeV

  float ee_rescaling = 1.035;
  float fh_rescaling = 1.095;
  float ah_rescaling = 1.095;
  
  char* outFileName = new char[1000];
  sprintf(outFileName,"./skimmed_ntuple_%s_needtoDelete.root", data);//,Min_, Max_); // the name of the file where you write the tree. 
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
      pi_event_orig = 0;
      pi_run = 0;
      pi_pdgID = 0;
      pi_beamEnergy = -1;
      pi_trueBeamEnergy = -1;
      pi_NRechits = 0;
      pi_NRechits_orig = 0;
      pi_combined_NRechits=0;
      pi_combined_rechits_energy.clear();
      pi_combined_rechit_x.clear();
      pi_combined_rechit_y.clear();
      pi_combined_rechit_z.clear();

      pi_rechit_energy.clear();
      pi_rechit_module.clear();
      pi_rechit_layer.clear();
      pi_rechit_chip.clear();
      pi_rechit_channel.clear();
      pi_rechit_type.clear();
      pi_rechit_x.clear();
      pi_rechit_y.clear();
      pi_rechit_z.clear();
      pi_rechit_iu.clear();
      pi_rechit_iv.clear();
      pi_rechit_iU.clear();
      pi_rechit_iV.clear();
      pi_rechit_amplitudeHigh.clear();
      pi_rechit_amplitudeLow.clear();
      pi_rechit_noise_flag.clear();
      pi_rechit_shower_start_layer = -999;
      pi_rechit_energyPerLayer=0;//.clear();
      pi_rechit_nHitsPerLayer=0;//.clear();
      pi_rechit_cogX=0;//.clear();
      pi_rechit_cogY=0;//.clear(); 
      // cout<<"insides the event loop: check1"<<endl;

      pi_ntracks = 0;
      pi_trackChi2_X = -1;
      pi_trackChi2_Y = -1;
      pi_dwcReferenceType = -1;
      pi_m_x = -999;
      pi_m_y = -999;
      pi_b_x = -999;
      pi_b_y = -999;
      pi_impactX_layer=0;//.clear();
      pi_impactY_layer=0;//.clear();
    
      pi_ahc_nHits = 0;
      pi_ahc_nHits_orig = 0;
      pi_ahc_energyPerLayer=0;// .clear();
      pi_ahc_nHitsPerLayer=0;//.clear();
      pi_ahc_hitI.clear();
      pi_ahc_hitJ.clear();
      pi_ahc_hitK.clear();
      pi_ahc_hitEnergy.clear();
      pi_ahc_hitX.clear();
      pi_ahc_hitY.clear();
      pi_ahc_hitZ.clear();
  
      pi_isHGC_AHC_sync = 0;
      pi_isGoodTrack = 0;
      pi_isFHNoisy = 0;
      pi_MuonVeto = 0;
      pi_isInTrackWindow = 0;
      pi_hgc_channel_mask=0;//->clear();
      pi_ahc_channel_mask=0;//->clear();
      pi_pass_noise_thres=0;//->clear();

      pi_energyLostEE = 0;
      pi_energyLostFH = 0;
      pi_energyLostBH = 0;
      pi_energyLostBeam = 0;
      pi_trimmed_ahcal_flag.clear();
      pi_rechitEn_trimAhcal.clear();
      pi_comb_rechit_x_trimAhcal.clear();
      pi_comb_rechit_y_trimAhcal.clear();
      pi_comb_rechit_z_trimAhcal.clear();
      pi_Nrechit_trimAhcal=0;;

      pi_lambda_trimAhcal.clear();
      pi_lambda.clear();
      event_count[0]++;
      pi_event_orig=event;
      pi_ahc_nHits_orig=ahc_nHits;
      pi_NRechits_orig=NRechits;
      pi_isHGC_AHC_sync=isHGC_AHC_sync;
      pi_isGoodTrack = isGoodTrack;
      pi_isFHNoisy = isFHNoisy;
      pi_MuonVeto=MuonVeto;
      pi_isInTrackWindow=isInTrackWindow;
      pi_hgc_channel_mask=hgc_channel_mask;
      pi_ahc_channel_mask=ahc_channel_mask;
      pi_pass_noise_thres=pass_noise_thres;
      
      event_count[1]++;
      if(!pi_isHGC_AHC_sync)
	{ cout<<"wrongevents"<<endl;}
      //if(Nrechits_FH_module_42 >80 ||  Nrechits_FH_module_45 > 80){count_fh++; continue;}
      if(isFHNoisy) {count_fh++; continue;}
      event_count[2]++;
      h_true_beamenergy[1]->Fill(trueBeamEnergy);
       //good track      
      bool isGoodTrack = (dwcReferenceType >= 13 && trackChi2_X < 10 && trackChi2_Y < 10);
      if(!isGoodTrack) { count_badtrack++; continue;}
      event_count[3]++;
      h_true_beamenergy[2]->Fill(trueBeamEnergy);
      //muon veto
      if(MuonVeto) continue;
      event_count[4]++;
      h_true_beamenergy[3]->Fill(trueBeamEnergy);
      //preshower
      if(rechit_shower_start_layer<=2) continue;
      event_count[5]++;
      h_true_beamenergy[4]->Fill(trueBeamEnergy);
      //track impact only for TB data & discrete simulation but not for flat energy samples (no track impact information was available)
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
      pi_rechit_cogX=rechit_cogX;
      pi_rechit_cogY=rechit_cogY;
      pi_impactX_layer=TrackImpactX_layer;
      pi_impactY_layer=TrackImpactY_layer;
      pi_rechit_energyPerLayer=rechit_energyPerLayer;
      pi_rechit_nHitsPerLayer=rechit_nHitsPerLayer;
      pi_rechit_shower_start_layer=rechit_shower_start_layer;
      pi_ntracks = ntracks;
      pi_trackChi2_X = trackChi2_X;
      pi_trackChi2_Y = trackChi2_Y;
      pi_dwcReferenceType = dwcReferenceType;
      pi_m_x = m_x;
      pi_m_y = m_y;
      pi_b_x = b_x;
      pi_b_y = b_y;

      //pi_ahc_nHits = ahc_nHits;
      pi_ahc_energyPerLayer = ahc_energyPerLayer;
      pi_ahc_nHitsPerLayer = ahc_nHitsPerLayer;
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
	  pi_lambda.push_back(lam_);
	  pi_lambda_trimAhcal.push_back(lam_);
	  float recx = rechit_x->at(i);
	  float	recy = rechit_y->at(i);
	  float recx_ini = rechit_x->at(i);
          float recy_ini = rechit_y->at(i);
	  float dX=0.0,dY=0.0,dX_in=0.0,dY_in=0.0;

	  //pi_trimmed_ahcal_flag.push_back(1);
	  if(!strcmp(data,"data")) {
	    //cout<<"inside the loop"<<endl;
	    if(layer<28 && layer==1)
              {
                h_xVs_trackX[0]->Fill(trackx,recx);//,trackx);
                h_yVs_trackY[0]->Fill(tracky,recy);
              }
            else
              {
		if(layer==29)
		  {
		    h_xVs_trackX[1]->Fill(trackx,recx);//,trackx);
		    h_yVs_trackY[1]->Fill(tracky,recy);//,tracky);
		  }
	      }
	    
	    trackx = -1*trackx;
	    tracky = -1*tracky;
	     if(layer<28 && layer==1)
              {
                h_xVs_1trackX[0]->Fill(trackx,recx);
                h_yVs_1trackY[0]->Fill(tracky,recy);
              }
            else
              {
		if(layer==29)
		  {
		    h_xVs_1trackX[1]->Fill(trackx,recx);
		    h_yVs_1trackY[1]->Fill(tracky,recy);
		  }
              }


	     dX = recx - trackx;
	     dY = recy - tracky;
	     dX_in= dX;
	     dY_in= dY;
	    if(layer==1)
	      {
		h_dXvsdY_EE[0]->Fill(dX,dY);
		h_dX_EE[0]->Fill(dX);
		h_dY_EE[0]->Fill(dY);
	      }
	    else if (layer==14)
	    {
	      h_dXvsdY_EE[1]->Fill(dX,dY);
	      h_dX_EE[1]->Fill(dX);
	      h_dY_EE[1]->Fill(dY);
	    }
          else if (layer==29)
	    {
	      h_dXvsdY_FH[0]->Fill(dX,dY);
	      h_dX_FH[0]->Fill(dX);
              h_dY_FH[0]->Fill(dY);
	    }
          else if(layer==34)
	    {
	      h_dXvsdY_FH[1]->Fill(dX,dY);
	      h_dX_FH[1]->Fill(dX);
              h_dY_FH[1]->Fill(dY);
	    }
	    h_muo_dxvsdy->Fill(dX,dY);
	    float dR = sqrt((dX*dX)+(dY*dY));
	    //h_dR->Fill(dR);
	    if(dr_a>dR) //to fill the minimum dR
	      dr_a=dR;
	    std::pair<float, float> dxy_al = dxy_alignment(layer);
	    float dx_corr = dxy_al.first;
	    float dy_corr = dxy_al.second; 
	    recx -= dx_corr; //correcting for alignment
	    recy -= dy_corr;
	    dX = recx - trackx;
	    dY = recy - tracky;
	    if(layer==1)	      
	      h_dXvsdY_EE_Acorr[0]->Fill(dX,dY);
	      
	    else if (layer==14)
	      h_dXvsdY_EE_Acorr[1]->Fill(dX,dY);
	    else if (layer==29)
	      h_dXvsdY_FH_Acorr[0]->Fill(dX,dY);
	    else if(layer==6)
	      h_dXvsdY_FH_Acorr[1]->Fill(dX,dY);

	    dR = sqrt((dX*dX)+(dY*dY));
	    if(dr_b>dR)
	      dr_b=dR;
	    h_muo_dxvsdy_afterCorr->Fill(dX,dY);
	    //h_dR_afterCorr->Fill(dR);

	    
	  }
	  float recx_fin = recx;
          float recy_fin = recy;
	  // regularizing the hits going beyond HGCAl physical boundaries
	  if(layer<=28)
	    {
	      if(recx_fin>7.5 || recx_fin <(-7.5) ||recy_fin>7.5 || recy_fin<(-7.5))
		{
		  count_rechit[1]++;
		  if(recx_fin>7.5)
		    recx = 7.5;
		  else if (recx_fin <(-7.5))
		    recx=-7.5;
		  else if(recy_fin>7.5)
		    recy=7.5;
		  else if(recy_fin<(-7.5))
		    recy=-7.5;
		      
		}
	    }
	  else
	    {
	      if(recx_fin>18.5 || recx_fin <(-18.5) ||recy_fin>18.5 || recy_fin<(-18.5))
		{
		  count_rechit[2]++;
		  if(recx_fin>18.5)
                    recx = 18.5;
                  else if (recx_fin <(-18.5))
                    recx=-18.5;
                  else if(recy_fin>18.5)
                    recy=18.5;
                  else if(recy_fin<(-18.5))
                    recy=-18.5;
		}

	    }

	  pi_rechit_layer.push_back(layer);
	  pi_rechit_energy.push_back(energy);
  	  pi_rechit_x.push_back(recx);//hit_x->at(i));
  	  pi_rechit_y.push_back(recy);//hit_y->at(i));
	  pi_rechitEn_trimAhcal.push_back(energy);
	  pi_comb_rechit_x_trimAhcal.push_back(recx);
	  pi_comb_rechit_y_trimAhcal.push_back(recy);
	  pi_comb_rechit_z_trimAhcal.push_back(rechit_z->at(i));
  	  pi_rechit_z.push_back(rechit_z->at(i));
  	  pi_rechit_iu.push_back(rechit_iu->at(i));
  	  pi_rechit_iv.push_back(rechit_iv->at(i));
  	  pi_rechit_iU.push_back(rechit_iU->at(i));
  	  pi_rechit_iV.push_back(rechit_iV->at(i));
  	  pi_rechit_amplitudeHigh.push_back(rechit_amplitudeHigh->at(i));
  	  pi_rechit_amplitudeLow.push_back(rechit_amplitudeLow->at(i));
  	  pi_rechit_noise_flag.push_back(rechit_noise_flag->at(i));
	  pi_rechit_module.push_back(rechit_module->at(i));
  	  pi_rechit_chip.push_back(rechit_chip->at(i));
  	  pi_rechit_channel.push_back(rechit_channel->at(i));
  	  pi_rechit_type.push_back(rechit_type->at(i));

	  pi_combined_rechits_energy.push_back(energy);
	  pi_combined_rechit_x.push_back(recx);//hit_x->at(i));
	  pi_combined_rechit_y.push_back(recy);//hit_y->at(i));
	  pi_combined_rechit_z.push_back(rechit_z->at(i));
	  nHgc++;
	  nrechits++;
	  nrechit_trimAhcal++;
	  h_rechitvslambda->Fill(lam_,energy);
	  //filling histograms energy wise
	  for(int ibin=0;ibin<8;ibin++)
	    {
	      float True_Energy =0.0;
	      if(!strcmp(data,"data"))
		True_Energy = beamEnergy;
	      else
		True_Energy =trueBeamEnergy;
	      if(True_Energy>=En_bin[ibin]-2 && True_Energy<=En_bin[ibin]+2)
		{		  
		  h_xVsy[ibin]->Fill(recx_ini,recy_ini);
		  h_xVsy_corr[ibin]->Fill(recx_fin,recy_fin);
		  if(layer<=28)
		    {
		      
		      h_xVsy_EE[ibin]->Fill(recx_ini,recy_ini);
		      h_xVsy_corr_EE[ibin]->Fill(recx_fin,recy_fin);
		      if(layer==1)
			{
			  h_xVsy_EE_l1[ibin]->Fill(recx_ini,recy_ini);
			  h_xVsy_EE_l1_Acorr[ibin]->Fill(recx_fin,recy_fin);
			}
		      else if(layer==14)
			{
			  h_xVsy_EE_l14[ibin]->Fill(recx_ini,recy_ini);
			 h_xVsy_EE_l14_Acorr[ibin]->Fill(recx_fin,recy_fin);
			}

		      if(rechit_shower_start_layer>20 && layer==1)
			{
			  h_xVsy_SS20_EE_l1[ibin]->Fill(dX_in,dY_in);//recx_ini,recy_ini);
			  h_xVsy_SS20_EE_l1_Acorr[ibin]->Fill(dX,dY);//recx_fin,recy_fin);
			}
		      else if (rechit_shower_start_layer>20 && layer==14)
			{
			  h_xVsy_SS20_EE_l14[ibin]->Fill(dX_in,dY_in);//recx_ini,recy_ini);
                          h_xVsy_SS20_EE_l14_Acorr[ibin]->Fill(dX,dY);//recx_fin,recy_fin);

			}
		    }
		  else
		    {
		      h_xVsy_FH[ibin]->Fill(recx_ini,recy_ini);
		      h_xVsy_corr_FH[ibin]->Fill(recx_fin,recy_fin);
		       if(layer==29)
                        {
                          h_xVsy_FH_l1[ibin]->Fill(recx_ini,recy_ini);
                          h_xVsy_FH_l1_Acorr[ibin]->Fill(recx_fin,recy_fin);
                        }
                      else if(layer==34)
                        {
                          h_xVsy_FH_l6[ibin]->Fill(recx_ini,recy_ini);
                         h_xVsy_FH_l6_Acorr[ibin]->Fill(recx_fin,recy_fin);
                        }

		    }

		}
	    }
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
      pi_NRechits = nHgc;

      
      // ////////////////////////////////////////////
      // //            AHCAL Part                  //
      // ////////////////////////////////////////////
      
      /// Read AHCAL Tree
      
      if(DEBUG) cout<<"DEBUG: Start Analylizing AHCAL  RecHits!!"<<endl;
      //cout<<ahc_nHits<<"ahcal"<<endl;
      float rechitEnergySum_AH_trimm=0.0;
      for(int i = 0 ; i < ahc_nHits; i++) {
      	if(ahc_channel_mask->at(i)) continue;
  	pi_ahc_hitI.push_back(ahc_hitI->at(i));
  	pi_ahc_hitJ.push_back(ahc_hitJ->at(i));
  	pi_ahc_hitK.push_back(ahc_hitK->at(i));
	float lam_ = lambda[ahc_hitK->at(i)+39];
	// float lam_pi = lambda_pi[ahc_hitK->at(i)+39];
	// pi_lambda_pi.push_back(lam_pi);
	pi_lambda.push_back(lam_);
	h_rechitvslambda->Fill(lam_,ahc_hitEnergy->at(i));
	count_rechit[0]++;
 	pi_ahc_hitZ.push_back(ahc_hitZ->at(i));
  	pi_ahc_hitEnergy.push_back(ahc_hitEnergy->at(i));
	int layer = ahc_hitK->at(i)+40;
	//alignment corrections
	double trackx = TrackImpactX_layer->at(layer-1);// track_x[layer-1];
	double tracky =  TrackImpactY_layer->at(layer-1);//track_y[layer-1];
	
	float recx =ahc_hitX->at(i);
	float recy =ahc_hitY->at(i);
	float recx_ini = recx;
	float recy_ini = recy;
	if(layer==41)
	  {
	    h_xVs_trackX[2]->Fill(trackx,recx);
	    h_yVs_trackY[2]->Fill(tracky,recy);
	  }
	//layer=ahc_hitK->at(i);
	if(!strcmp(data,"data")) {
	  // trackx = -1*trackx;
	  // tracky = -1*tracky;
	  if(layer==41)
	    {
	      h_xVs_1trackX[2]->Fill(trackx,recx);
	      h_yVs_1trackY[2]->Fill(tracky,recy);
	    }
	  float dX = recx - trackx;
	  float dY = recy - tracky;
	  
	  if(layer-40==1)
	    {
	    h_dXvsdY_AH[0]->Fill(dX,dY);
	    h_dX_AH[0]->Fill(dX);
	    h_dY_AH[0]->Fill(dY);
	    }
	  else if (layer-40==20)
	    {
	      h_dXvsdY_AH[1]->Fill(dX,dY);
	      h_dX_AH[1]->Fill(dX);
	      h_dY_AH[1]->Fill(dY);

	    }

	  h_muo_dxvsdy->Fill(dX,dY);
	  float dR = sqrt((dX*dX)+(dY*dY));
	  if(dr_a>dR)
              dr_a=dR;

	  //h_dR->Fill(dR);
	  std::pair<float, float> dxy_al = dxy_alignment(layer);
	  float dx_corr = dxy_al.first;
	  float dy_corr = dxy_al.second;
	  recx -= dx_corr;
	  recy -= dy_corr;
	  dX = recx - trackx;
	  dY = recy - tracky;
	   if(layer-40==1)
            h_dXvsdY_AH_Acorr[0]->Fill(dX,dY);
          else if (layer-40==20)
            h_dXvsdY_AH_Acorr[1]->Fill(dX,dY);

	  dR = sqrt((dX*dX)+(dY*dY));
	  if(dr_b>dR)
              dr_b=dR;

	  h_muo_dxvsdy_afterCorr->Fill(dX,dY);
	  //h_dR_afterCorr->Fill(dR);

	  // std::pair<float, float> dxy_al = dxy_alignment(layer);
	  // float dx_corr = dxy_al.first;
	  // float dy_corr = dxy_al.second;
	  // recx -= dx_corr;
	  // recy -= dy_corr;
	  }
	float recx_fin = recx;
        float recy_fin = recy;
	//regularizning the alignment
       if(recx_fin>36 || recx_fin <(-36) ||recy_fin>36 || recy_fin<(-36))
	 {
                count_rechit[3]++;
		if(recx_fin>18.5)
                    recx = 18.5;
                  else if (recx_fin <(-18.5))
                    recx=-18.5;
                  else if(recy_fin>18.5)
                    recy=18.5;
                  else if(recy_fin<(-18.5))
                    recy=-18.5;


	 }
	for(int ibin=0;ibin<8;ibin++)
            {
	      float True_Energy =0.0;
              if(!strcmp(data,"data"))
                True_Energy = beamEnergy;
              else
                True_Energy =trueBeamEnergy;
              if(True_Energy>=En_bin[ibin]-2 && True_Energy<=En_bin[ibin]+2)
                {
                  h_xVsy[ibin]->Fill(recx_ini,recy_ini);
                  h_xVsy_corr[ibin]->Fill(recx_fin,recy_fin);
		  h_xVsy_AH[ibin]->Fill(recx_ini,recy_ini);
		  h_xVsy_corr_AH[ibin]->Fill(recx_fin,recy_fin);
		   if(layer-40==1)
			{
                          h_xVsy_AH_l1[ibin]->Fill(recx_ini,recy_ini);
                          h_xVsy_AH_l1_Acorr[ibin]->Fill(recx_fin,recy_fin);
                        }
                      else if(layer-40==20)
                        {
                          h_xVsy_AH_l20[ibin]->Fill(recx_ini,recy_ini);
                         h_xVsy_AH_l20_Acorr[ibin]->Fill(recx_fin,recy_fin);
                        }

                }
            }

	h_dR_afterCorr->Fill(dr_b);
	h_dR->Fill(dr_a);
	pi_ahc_hitX.push_back(recx);
        pi_ahc_hitY.push_back(recy);

	pi_combined_rechits_energy.push_back(ahc_hitEnergy->at(i));
	pi_combined_rechit_x.push_back(recx);
	pi_combined_rechit_y.push_back(recy);
	pi_combined_rechit_z.push_back(ahc_hitZ->at(i)+offset);
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
		pi_trimmed_ahcal_flag.push_back(1);
		//pi_lambda_pi_trimAhcal.push_back(lam_pi);
		pi_lambda_trimAhcal.push_back(lam_);
		pi_rechitEn_trimAhcal.push_back(ahc_hitEnergy->at(i));
		pi_comb_rechit_x_trimAhcal.push_back(recx);
		pi_comb_rechit_y_trimAhcal.push_back(recy);
		pi_comb_rechit_z_trimAhcal.push_back(ahc_hitZ->at(i)+offset);
		rechitEnergySum_AH_trimm+=ahc_hitEnergy->at(i);
		nrechit_trimAhcal++;
		count++;		
	      }
	  }
	//cout<<count<<endl;
	//if(count==0)
	//	pi_trimmed_ahcal_flag.push_back(count);
	  
      }// AHCAL hit loop
   
      // code to get the updated relative weight between FH and AH after emoving additional layers
      if(trueBeamEnergy>=48 && trueBeamEnergy<=52)
        {
	  if(rechit_shower_start_layer>28)
	    {
	      float rel_weight = 0.02;
	      float step_ = rel_weight;
	      for(int ii = 0; ii < 50; ii++)
		{
		  rel_weight = step_*(ii+1);
		  h_rechit_energy_FB_rel_weightScan[ii]->Fill(Esum_rechits_FH+(rel_weight*rechitEnergySum_AH_trimm)); 
		}
	      h_mipsinEE->Fill((Esum_rechits_EE*0.0105)+(0.0789*Esum_rechits_FH)+(0.07574*rechitEnergySum_AH_trimm));
	    }
	}
      h_true_beamenergy[5]->Fill(trueBeamEnergy);
      pi_ahc_nHits=nAhc;
      pi_combined_NRechits=nrechits;
      pi_Nrechit_trimAhcal=nrechit_trimAhcal;
      nrechit_trimAhcal=0;
      pion_tree->Fill();
      // cout<<"total rechits"<<"\t"<<count_rechit[0]<<"\t"<<count_rechit[1]<<"\t"<<count_rechit[2]<<endl;
      // cout<<"fraction"<<"\t"<<count_rechit[1]/count_rechit[0]<<"\t"<<count_rechit[2]/count_rechit[0]<<endl;
      hist_total->Fill(count_rechit[0]);
      float frac = count_rechit[1];///count_rechit[0];
      hist_frac_EE->Fill(frac);
      frac = count_rechit[2];///count_rechit[0];
      hist_frac_FH->Fill(frac);
      frac = count_rechit[3];///count_rechit[0];
      hist_frac_AH->Fill(frac);
      
      count_rechit[0]=0;
      count_rechit[1]=0;
      count_rechit[2]=0;
      count_rechit[3]=0;
      rechits_AH = nAhc;
      total_rechits = nrechits;
      nrechits=0;
      nHgc=0;
      nAhc =0;
      Esum_rechits_EE_inGeV=0.0105*Esum_rechits_EE;
      Esum_rechits_FH_inGeV= 0.0789*Esum_rechits_FH;
      Esum_rechits_AH_inGeV=0.0316*Esum_rechits_AH;
      if((trueBeamEnergy/(Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV))>4) //continue;
	{
      h_nrechits->Fill(total_rechits);
      h_nrechits_EE->Fill(rechits_EE);
      h_nrechits_FH->Fill(rechits_FH);
      h_nrechits_AH->Fill(rechits_AH);
      h_nrechits_vs_EE->Fill(total_rechits,rechits_EE);
      h_nrechits_vs_FH->Fill(total_rechits,rechits_FH);
      h_nrechits_vs_AH->Fill(total_rechits,rechits_AH);
      h_nrechits_EE_vs_FH->Fill(rechits_EE, rechits_FH);
      nAhc =0;
      h_EnergySum_inEE->Fill(Esum_rechits_EE);
      h_EnergySum_inFH->Fill( Esum_rechits_FH);
      h_EnergySum_inEE_inGeV->Fill(0.0105*1.5*Esum_rechits_EE);
      h_EnergySum_inFH_inGeV->Fill(0.0789*Esum_rechits_FH);
      h_EnergySum_inEE_Mipsvs_true->Fill(Esum_rechits_EE,trueBeamEnergy);
      h_EnergySum_inFH_Mipsvs_true->Fill(Esum_rechits_FH,trueBeamEnergy);

      // nrechits=0;
      // nHgc=0;
      // Esum_rechits_EE_inGeV=0.0105*Esum_rechits_EE;
      // Esum_rechits_FH_inGeV= 0.0789*Esum_rechits_FH;
      // Esum_rechits_AH_inGeV=0.0316*Esum_rechits_AH;
      // if(trueBeamEnergy>=48 && trueBeamEnergy<=52)
      // 	{
      // if(rechit_shower_start_layer<=28)
      // 	h_ssinEE->Fill(Esum_rechits_EE_inGeV + Esum_rechits_FH_inGeV);//+Esum_rechits_AH_inGeV);
      // else//(rechit_shower_start_layer>28)
      // 	h_mipsinEE->Fill(Esum_rechits_EE_inGeV + Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV);
      // 	}
      h_EnergySum_Rechit->Fill(Esum_rechits_EE+Esum_rechits_FH+Esum_rechits_AH);
      h_EnergySum_inGeV->Fill((Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV));
      h_EnergySum_inAH->Fill(Esum_rechits_AH);
      h_EnergySum_inAH_inGeV->Fill(Esum_rechits_AH_inGeV);
      h_EnergySum_inGeVvstrue->Fill(Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV, trueBeamEnergy);
      h_EnergySum_inMipsvstrue->Fill(Esum_rechits_EE+Esum_rechits_FH+Esum_rechits_AH, trueBeamEnergy);
      h_EnergySum_inAH_Mipsvs_true->Fill(Esum_rechits_AH,trueBeamEnergy);
      h_EnergySum_inAH_vs_true->Fill(Esum_rechits_AH_inGeV,trueBeamEnergy);
      h_EnergySum_inFH_vs_true->Fill(Esum_rechits_FH_inGeV,trueBeamEnergy);
      h_EnergySum_inEE_vs_true->Fill(Esum_rechits_EE_inGeV,trueBeamEnergy);
      h_EnergySum_ratio_inGeV->Fill(trueBeamEnergy/(Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV));
      h_EnergySum_ratio_flipped->Fill((Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV)/trueBeamEnergy);
      h_EnergySum_ratio_inMips->Fill(trueBeamEnergy/(Esum_rechits_EE+Esum_rechits_FH+Esum_rechits_AH));
      h_EnergySum_inEE_vs_FH->Fill(Esum_rechits_EE,Esum_rechits_FH);
      h_EnergySum_inEE_vs_FH_inGeV->Fill(Esum_rechits_EE_inGeV,Esum_rechits_FH_inGeV);
      h_EnergySum_ratio_vs_true->Fill(trueBeamEnergy/(Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV), trueBeamEnergy);
      h_EnergySum_ratio_vs_true_flipped->Fill((Esum_rechits_EE_inGeV+Esum_rechits_FH_inGeV+Esum_rechits_AH_inGeV)/trueBeamEnergy, trueBeamEnergy);
      }

      double rechitEnergySum_EE = Esum_rechits_EE;
      double rechitEnergySum_FH = Esum_rechits_FH;
      double rechitEnergySum_AH = Esum_rechits_AH;
      double EE_detscale = Esum_rechits_EE_inGeV;//(0.0105rechitEnergySum_EE/EE_scale);                                        
      double FH_detscale = Esum_rechits_FH_inGeV;//(rechitEnergySum_FH/FH_AH_scale);                                                 
      double AH_detscale = Esum_rechits_AH_inGeV;//(alpha_*rechitEnergySum_AH)/FH_AH_scale;                                   
      double full_energy = FH_detscale+AH_detscale;                                                                                       
      double total_energy= (rechitEnergySum_EE/EE_scale) + (rechitEnergySum_FH+ (alpha_ * rechitEnergySum_AH))/FH_AH_scale;
      double tot_E_gev= total_energy;
      ///////////////////////////////////////////////////////////                                                                           
      /////     H hadrons ; Chi2 matrix initialzation     /////                                                                             
      //////////////////////////////////////////////////////////
      //filling the chi2 energy
      //cout<<"check2"<<"H-hadronsinitialize tree branch"<<endl;
      total_energy = (rechitEnergySum_FH+ (alpha_ * rechitEnergySum_AH))/FH_AH_scale+0.4;

      double E_beam =0.0;
      double O = 0.0;
      double sigma2_H = (0.089*0.089)*(total_energy*total_energy) + (1.25*1.25)*total_energy;
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
		  
		  //cout<<"check8"<<"in the loop"<<"\t"<<i_bin<<endl;
		  
		  float w1 = getChi2Weights_H(E_beam).at(0);
		  float w2 = getChi2Weights_H(E_beam).at(1);
		  float w3 = getChi2Weights_H(E_beam).at(2);
		  //cout<<"check5"<<"in the loop"<<"\t"<<i_bin<<endl;
		  float chi2_energy = w1*rechitEnergySum_EE + w2*rechitEnergySum_FH + w3*rechitEnergySum_AH + O;
		  if(chi2_method > 0)  {
		    if(chi2_method == 1) O = 0.4;
		    chi2_energy = w1*1.0 + w2*FH_detscale + w3*AH_detscale + O;
		  }
		  hist_resp_total_1[i_bin]->Fill(chi2_energy);
		  hist_resp_SS_FH_1[i_bin]->Fill(chi2_energy);
		  hist_fixwt_total[i_bin]->Fill(total_energy);
                  hist_fixwt_SS_FH[i_bin]->Fill(total_energy);

		  ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
		  ROOT::Math::SVector<double, 2> consts_H;
		  ROOT::Math::SVector<double, 2> values_H;
		  for(int i = 0; i < 2; i++) {
		    for(int j = 0; j < 2; j++) {
		      coeffs_H(i,j) = 0.0;
		    }	
		    consts_H(i) = 0.0;
		    values_H(i) = 0.0;
		  }
		  IsIn_en_range = (total_energy > en_range_FH_xmin_vec[i_bin][Elist[i_bin]+2] && total_energy < en_range_FH_xmax_vec[i_bin][Elist[i_bin]+2]);
		  //IsIn_en_range = true;
		  if(IsIn_en_range){

		  if(chi2_method == 0)
		    {
		      coeffs_H(0,0) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2_H;
		      coeffs_H(0,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2_H;
		      coeffs_H(1,0) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2_H;
		      coeffs_H(1,1) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2_H;
		      
		      consts_H(0)   = (E_beam-O)*rechitEnergySum_FH/sigma2_H;
		      consts_H(1)   = (E_beam-O)*rechitEnergySum_AH/sigma2_H;
		      consts_H_vec[i_bin]+= consts_H;
		      coeff_H_vec[i_bin]+=coeffs_H;
		    }
		  else
		    {
		      coeffs_H(0,0) = (FH_detscale*FH_detscale)/sigma2_H;
		      coeffs_H(0,1) = (AH_detscale*FH_detscale)/sigma2_H;
		      coeffs_H(1,0) = (FH_detscale*AH_detscale)/sigma2_H;
		      coeffs_H(1,1) = (AH_detscale*AH_detscale)/sigma2_H;
		      
		      consts_H(0)   = (E_beam-O)*FH_detscale/sigma2_H;
		      consts_H(1)   = (E_beam-O)*AH_detscale/sigma2_H;
		      consts_H_vec[i_bin] += consts_H;
		      coeff_H_vec[i_bin]+=coeffs_H;

		    }
		  }
		}
	      else if (rechit_shower_start_layer<=28)
		{
		  //cout<<"check6"<<"in the loop"<<"\t"<<i_bin<<endl;
		  float w1 = getChi2Weights_EH(E_beam).at(0);                                                                          
		  float w2 = getChi2Weights_EH(E_beam).at(1);                                                                         
		  float w3 = getChi2Weights_EH(E_beam).at(2);
		  //cout<<"check6"<<"in the loop"<<"\t"<<i_bin<<endl;
		  float chi2_energy = w1*rechitEnergySum_EE + w2*rechitEnergySum_FH + w3*rechitEnergySum_AH + O;
		  if(chi2_method > 0)  {
                    if(chi2_method == 1) O = 0.4;
                    chi2_energy = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale + O;
                  }

		  hist_resp_total_1[i_bin]->Fill(chi2_energy);
                  hist_resp_SS_EE_1[i_bin]->Fill(chi2_energy);
		  hist_fixwt_total[i_bin]->Fill(tot_E_gev);
                  hist_fixwt_SS_EE[i_bin]->Fill(tot_E_gev);


		  double sigma2_EH =(0.084*0.084)*(tot_E_gev*tot_E_gev) + (1.39*1.39)*tot_E_gev;// (0.084*0.084) + (1.39*1.39)*tot_E_gev;
		  ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
		  ROOT::Math::SVector<double, 3> consts_EH;
		  ROOT::Math::SVector<double, 3> values_EH;

		  for(int i = 0; i < 3; i++) {
		    for(int j = 0; j < 3; j++) {
		      coeffs_EH(i,j) = 0.0;
		    }
		    consts_EH(i) = 0.0;
		    values_EH(i) = 0.0;
		  }
		  
		   IsIn_en_range = (tot_E_gev > en_range_EE_xmin_vec[i_bin][Elist[i_bin]+2] && tot_E_gev < en_range_EE_xmax_vec[i_bin][Elist[i_bin]+2]);
		   //IsIn_en_range = true;
		   if(IsIn_en_range){
		     if(chi2_method==0)
		       {
		      coeffs_EH(0,0) = (rechitEnergySum_EE*rechitEnergySum_EE)/sigma2_EH;
		      coeffs_EH(0,1) = (rechitEnergySum_EE*rechitEnergySum_FH)/sigma2_EH;
		      coeffs_EH(0,2) = (rechitEnergySum_EE*rechitEnergySum_AH)/sigma2_EH;
		      coeffs_EH(1,0) = (rechitEnergySum_FH*rechitEnergySum_EE)/sigma2_EH;
		      coeffs_EH(1,1) = (rechitEnergySum_FH*rechitEnergySum_FH)/sigma2_EH;
		      coeffs_EH(1,2) = (rechitEnergySum_FH*rechitEnergySum_AH)/sigma2_EH;
		      coeffs_EH(2,0) = (rechitEnergySum_AH*rechitEnergySum_EE)/sigma2_EH;
		      coeffs_EH(2,1) = (rechitEnergySum_AH*rechitEnergySum_FH)/sigma2_EH;
		      coeffs_EH(2,2) = (rechitEnergySum_AH*rechitEnergySum_AH)/sigma2_EH;

		      consts_EH(0)   = E_beam*rechitEnergySum_EE/sigma2_EH;
		      consts_EH(1)   = E_beam*rechitEnergySum_FH/sigma2_EH;
		      consts_EH(2)   = E_beam*rechitEnergySum_AH/sigma2_EH;
		      consts_EH_vec[i_bin]+= consts_EH;
		      coeff_EH_vec[i_bin]+=coeffs_EH;

		    }
		  else
		    {
		      coeffs_EH(0,0) = (EE_detscale*EE_detscale)/sigma2_EH;
		      coeffs_EH(0,1) = (EE_detscale*FH_detscale)/sigma2_EH;
		      coeffs_EH(0,2) = (EE_detscale*AH_detscale)/sigma2_EH;
		      coeffs_EH(1,0) = (FH_detscale*EE_detscale)/sigma2_EH;
		      coeffs_EH(1,1) = (FH_detscale*FH_detscale)/sigma2_EH;
		      coeffs_EH(1,2) = (FH_detscale*AH_detscale)/sigma2_EH;
		      coeffs_EH(2,0) = (AH_detscale*EE_detscale)/sigma2_EH;
		      coeffs_EH(2,1) = (AH_detscale*FH_detscale)/sigma2_EH;
		      coeffs_EH(2,2) = (AH_detscale*AH_detscale)/sigma2_EH;


		      consts_EH(0)   = (E_beam - O)*EE_detscale/sigma2_EH;
		      consts_EH(1)   = (E_beam - O)*FH_detscale/sigma2_EH;
		      consts_EH(2)   = (E_beam - O)*AH_detscale/sigma2_EH;
		      consts_EH_vec[i_bin]+= consts_EH;
		      coeff_EH_vec[i_bin]+=coeffs_EH;
		      
		    }
		   }
		}
            }	         
        }//true bin loop
      //cout<<"check3"<<"After chi2"<<endl;
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      
      //if(DEBUG && jentry > 10) break;
      //      gSystem->Exit(0);
      //if(jentry > 2) break;
      //pion_tree->Fill();
    } // loop over entries
  
  char* name = new char[1000];
  sprintf(name,"chi2_calibFact_H_hadrons_flatEn_%d_scalMC_2sigma.txt",chi2_method);
  std::ofstream file_H;
  file_H.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_H_hadrons_flatEn_%d.txt",chi2_method,ios::out);                        
  sprintf(name,"chi2_calibFact_EH_hadrons_flatEn_%d_scalMC_2sigma.txt",chi2_method);
  std::ofstream file_EH;
  file_EH.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_EH_hadrons_flatEn_%d.txt",chi2_method,ios::out

  /////////////////////////////////////////////////////                                                                                     
  //////                 EH Hadrons            ////////                                                                                     
  /////////////////////////////////////////////////////                                                                                     

   for(int i_en =0; i_en<85; i_en++)
    {
      ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
      ROOT::Math::SVector<double, 3> consts_EH;
      ROOT::Math::SVector<double, 3> values_EH;
      ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
      ROOT::Math::SVector<double, 2> consts_H;
      ROOT::Math::SVector<double, 2> values_H;
       //      bool isInverted_EH = false;                                                                                                  
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          coeffs_EH(i,j) = 0.0;
        }
        consts_EH(i) = 0.0;
        values_EH(i) = 0.0;
      }
       for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
          coeffs_H(i,j) = 0.0;
        }
        consts_H(i) = 0.0;
        values_H(i) = 0.0;
      }

      consts_EH = consts_EH_vec[i_en];
      values_EH = values_EH_vec[i_en];
      coeffs_EH = coeff_EH_vec[i_en];
      isInverted_EH = coeffs_EH.Invert();
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          cout<<coeffs_EH(i, j)<<"\t";
        }
        cout<<endl;
      }

      cout<<" "<<consts_EH(0)<<"\t"<<consts_EH(1)<<"\t"<<consts_EH(2)<<endl;
consts_H = consts_H_vec[i_en];
      values_H = values_H_vec[i_en];
      coeffs_H = coeff_H_vec[i_en];

      isInverted_H = coeffs_H.Invert();
      for(int i = 0; i <3; i++) {
        for(int j = 0; j < 3; j++) {
          cout<<coeffs_H(i, j)<<"\t";
        }
        cout<<endl;
      }
      cout<<" "<<consts_H(0)<<"\t"<<consts_H(1)<<"\t"<<consts_H(2)<<endl;
      if(isInverted_EH) {
        values_EH = coeffs_EH*consts_EH;
        cout<<"EH Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = "<<values_EH(0)<<" ;w2 = "<<values_EH(1)<<" ;w3 = "<<values_EH(2)<<";"<<endl;
	cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_EH(0,0))<< " " << sqrt(coeffs_EH(1,1)) << " " << sqrt(coeffs_EH(2,2))<<endl;

         file_EH<<Elist[i_en]+2<<"\t"<<values_EH(0)<<"\t"<<values_EH(1)<<"\t"<<values_EH(2)<<"\t"<<sqrt(coeffs_EH(0,0))<<"\t"<<sqrt(coeffs_EH(1,1))<<"\t"<<sqrt(coeffs_EH(2,2))<<"\n";

        //file_EH<<Elist[i_en]+2<<"\t"<<values_EH(0)<<"\t"<<values_EH(1)<<"\t"<<values_EH(2)<<"\n";
      }
      else {
        cout<<"Error: Could not invert for EH hadrons..."<<endl;
      }
      //isInverted_H = true;                                                                                                                
      if(isInverted_H) {
        values_H = coeffs_H*consts_H;
        //cout<<"H Hadrons => For E = "<<E_beam<<"GeV, w1 = "<<values_H(0)<<" ;w2 = "<<values_H(1)<<" ;w3 = "<<values_H(2)<<";"<<endl;      
        cout<<"H Hadrons => For E = "<<Elist[i_en]+2<<"GeV, w1 = 0.0 ;w2 = "<<values_H(0)<<" ;w3 = "<<values_H(1)<<";"<<endl;
	cout<<">>>>>>> Uncertainity  : "<<sqrt(coeffs_H(0,0)) << " " <<sqrt(coeffs_H(1,1))<<endl;

        file_H<<Elist[i_en]+2<<"\t"<<0<<"\t"<<values_H(0)<<"\t"<<values_H(1)<<"\t"<<0<<"\t"<<sqrt(coeffs_H(0,0))<<"\t"<<sqrt(coeffs_H(1,1))<<"\n";

        //file_H<<Elist[i_en]+2<<"\t"<<values_H(0)<<"\t"<<values_H(1)<<"\n";
      }
      else {
        cout<<"Error: Could not invert for H hadrons..."<<endl;
      }
    }


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
