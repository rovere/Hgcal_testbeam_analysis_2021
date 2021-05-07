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




int main(int argc, char* argv[])
{

  if (argc < 6) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" " << "configuration" 
	 <<" " << "energy" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *config          = argv[4];
  const char *energy = argv[5];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy);
  cout << "dataset " << data << " " << endl;
  cout << "configuration " << config << " " << endl;
  cout << "energy " << energy << " " << endl;

  hgcOctTB.EventLoop(data,energy);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data, const char *energy) {


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
  int nHgc=0, nAhc=0;
  //cout<<energy<<endl;
  char* outFileName = new char[1000];
  sprintf(outFileName,"./skimmed_ntuple_sim_config22_pion%sGeV_combinedHgc_Ahc_v46_patchMIP.root",energy);
  TFile* outfile = TFile::Open(outFileName,"recreate");
  //TTree* pion_tree;
  pion_tree = new TTree("pion_variables_v1","variables for pion analysis v1");
  HGCNtupleVariables::init_piTree();

  
  if(DEBUG) cout<<"DEBUG: Configuration = "<<conf_<<endl;
  if(DEBUG) cout<<"DEBUG: TOTAL_ACTIVE_LAYER = "<<TOTAL_ACTIVE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: EE_LAYER = "<<EE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: FH_LAYER = "<<FH_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: AH_LAYER = "<<AH_LAYER<<endl;

  if(DEBUG) cout<<"DEBUG: Entering event Loop"<<endl;

  Long64_t jentry;

  for (jentry=0; jentry<nentries;jentry++,hgc_jentry++)
  {
  //  for (jentry=0; jentry<1;jentry++,hgc_jentry++) {
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

      //cout<<"insides the event loop: check2"<<endl;

      //    pion_tree->Fill();
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
      //if(Nrechits_FH_module_42 >80 ||  Nrechits_FH_module_45 > 80){count_fh++; continue;}
      if(isFHNoisy) {count_fh++; continue;}
      event_count[2]++;
       //good track
      bool isGoodTrack = (dwcReferenceType >= 13 && trackChi2_X < 10 && trackChi2_Y < 10);
      if(!isGoodTrack) { count_badtrack++; continue;}
      event_count[3]++;
      
      //muon veto
      if(MuonVeto) continue;
      event_count[4]++;
      if(rechit_shower_start_layer<=2) continue;
      event_count[5]++;
      if(!isInTrackWindow) continue;
      event_count[6]++;
      h_beamenergy->Fill(beamEnergy);
      h_true_beamenergy->Fill(trueBeamEnergy);
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
      
      //pion_tree->Fill();

      
      ////////////////////////////////////////////
      //            HGCAL Part                  //
      ////////////////////////////////////////////
      
      /// Read HGCAL Tree
      if(DEBUG) cout<<"DEBUG: Start Analylizing HGCAL  RecHits!!"<<endl;
      // total_rechits_energy_EE=0;
      // total_rechits_energy_FH=0;
      for( int i = 0 ; i < NRechits; i++)
      	{
      	  int channel = rechit_channel->at(i);
      	  int chip = rechit_chip->at(i);
      	  int en_chan = 1000*chip + channel;
      	  int layer = rechit_layer->at(i);
      	  float energy=rechit_energy->at(i);
	  
  	  if(hgc_channel_mask->at(i)) continue;
  	  if(!pass_noise_thres->at(i))continue;
	  //cout<<rechit_x->at(i)<<"x"<<rechit_y->at(i)<<"y"<<endl;
	  pi_rechit_energy.push_back(energy);
  	  pi_rechit_x.push_back(rechit_x->at(i));
  	  pi_rechit_y.push_back(rechit_y->at(i));
  	  pi_rechit_z.push_back(rechit_z->at(i));
  	  pi_rechit_iu.push_back(rechit_iu->at(i));
  	  pi_rechit_iv.push_back(rechit_iv->at(i));
  	  pi_rechit_iU.push_back(rechit_iU->at(i));
  	  pi_rechit_iV.push_back(rechit_iV->at(i));
  	  pi_rechit_amplitudeHigh.push_back(rechit_amplitudeHigh->at(i));
  	  pi_rechit_amplitudeLow.push_back(rechit_amplitudeLow->at(i));
  	  pi_rechit_noise_flag.push_back(rechit_noise_flag->at(i));
	  pi_rechit_module.push_back(rechit_module->at(i));
  	  //pi_rechit_layer->push_back(rechit_layer->at(i));
  	  pi_rechit_chip.push_back(rechit_chip->at(i));
  	  pi_rechit_channel.push_back(rechit_channel->at(i));
  	  pi_rechit_type.push_back(rechit_type->at(i));
	  nHgc++;
  	  } //Nrechits loop
      //pion_tree->Fill();
      pi_NRechits = nHgc;
      
      // pi_rechit_x->clear();
      // pi_rechit_y->clear();
      // //////////////////////////////////////////
      
      
      // //    cout<<NRechits<<"\t"<<count<<endl;
      
      // ////////////////////////////////////////////
      // //            AHCAL Part                  //
      // ////////////////////////////////////////////
      
      /// Read AHCAL Tree
      if(DEBUG) cout<<"DEBUG: Start Analylizing AHCAL  RecHits!!"<<endl;

      for(int i = 0 ; i < ahc_nHits; i++) {
      	if(ahc_channel_mask->at(i)) continue;
  	pi_ahc_hitI.push_back(ahc_hitI->at(i));
  	pi_ahc_hitJ.push_back(ahc_hitJ->at(i));
  	pi_ahc_hitK.push_back(ahc_hitK->at(i));
  	pi_ahc_hitX.push_back(ahc_hitX->at(i));
  	pi_ahc_hitY.push_back(ahc_hitY->at(i));
  	pi_ahc_hitZ.push_back(ahc_hitZ->at(i));
  	pi_ahc_hitEnergy.push_back(ahc_hitEnergy->at(i));
  	// pi_ahc_energyPerLayer.push_back(ahc_energyPerLayer->at(i));
  	// pi_ahc_nHitsPerLayer.push_back(ahc_nHitsPerLayer->at(i));
	nAhc++;
      }//
      pi_ahc_nHits=nAhc;
      pion_tree->Fill();
      nAhc =0;
      nHgc=0;
      // h_total_rechitsEenrgy_nocut->Fill(total_rechits_energy_nocut);
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      
      //if(DEBUG && jentry > 10) break;
      //      gSystem->Exit(0);
      //      if(jentry > 10) break;
      //pion_tree->Fill();
    } // loop over entries
  

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
