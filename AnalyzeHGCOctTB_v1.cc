#define AnalyzeHGCOctTB_cxx
#include <iostream>
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>

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

  hgcOctTB.EventLoop(data);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data) {


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nentries3 = fChain3->GetEntriesFast();
  Long64_t hgc_jentry = 0;

  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  Long64_t nbytes3 = 0, nb3 = 0;

  

  int decade = 0;
  int ahc_zeroHit = 0;



  bool DEBUG = false;

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





  if(DEBUG) cout<<"DEBUG: Configuration = "<<conf_<<endl;
  if(DEBUG) cout<<"DEBUG: TOTAL_ACTIVE_LAYER = "<<TOTAL_ACTIVE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: EE_LAYER = "<<EE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: FH_LAYER = "<<FH_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: AH_LAYER = "<<AH_LAYER<<endl;

  if(DEBUG) cout<<"DEBUG: Entering event Loop"<<endl;

  Long64_t jentry;

  for (jentry=0; jentry<nentries;jentry++,hgc_jentry++)
    {
      //for (jentry=0; jentry<1;jentry++,hgc_jentry++) {
      // ==============print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;
    
      // ===============read this entry == == == == == == == == == == ==

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { break; cout<<"Breaking"<<endl;}
      nb = fChain->GetEntry(hgc_jentry);   nbytes += nb;
      nb2 = fChain2->GetEntry(hgc_jentry); nbytes2 += nb2;
      nb3 = fChain3->GetEntry(jentry); nbytes3 += nb3;

    
      // ===============Synch HGCAL & AHCAL with same entry == == == == == == == == == == ==

      while(run != runNumber)
	{
	  if(nentries>nentries3)
	    {
	
	      hgc_jentry++;
	      nb = fChain->GetEntry(hgc_jentry);   nbytes += nb;
	      nb2 = fChain2->GetEntry(hgc_jentry); nbytes2 += nb2;
	      cout<<"hgcal extra event"<<endl;
	    }
	  else if(nentries<nentries3)
	    {
	      jentry++;
	      nb3 = fChain3->GetEntry(jentry); nbytes3 += nb3;
	      cout<<"Ahcal extra event"<<endl;
	    }
	}


      h_beamenergy->Fill(beamEnergy);
      h_particle->Fill(pdgID);



      ////////////////////////////////////////////
      //            HGCAL Part                  //
      ////////////////////////////////////////////
      
      /// Read HGCAL Tree
      
      if(DEBUG) cout<<"DEBUG: Start Analylizing HGCAL  RecHits!!"<<endl;
      
      for(int i = 0 ; i < NRechits; i++)
	{
	  
	  if((*rechit_noise_flag)[i] == 0)
	    {
	      if((*rechit_chip)[i] !=3 && ((*rechit_channel)[i]!=28 || (*rechit_channel)[i]!=22) && (*rechit_chip)[i] !=4 && (*rechit_channel)[i]!=44 )
		{
	    
		  //cout<<(*rechit_chip)[i]<<"\t"<<(*rechit_channel)[i]<<"\t"<<(*rechit_type)[i]<<endl;
	    
		  if((*rechit_energy)[i]>0.5)
		    {
		      if((*rechit_layer)[i] ==1 && (*rechit_chip)[i] ==0)
			{
			  continue;
			}

		      else
			{
			  count_afterCuts++;
			  //cout<<(*rechit_detid)[i]<<"\t"<<(*rechit_module)[i]<<"\t"<<(*rechit_chip)[i]<<"\t"<<(*rechit_channel)[i]<<"\t"<<(*rechit_type)[i]<<endl;
			
			  rechits_layer[(*rechit_layer)[i]]+=(*rechit_energy)[i];
			  //bool alpha1 = (*rechit_noise_flag)[i];
			  //cout<<(*rechit_noise_flag)[i]<<endl;
			  
			  if((*rechit_layer)[i]<=28)
			    {
			      Esum_rechits_EE+=(*rechit_energy)[i];
			    }
			  else
			    {
			      Esum_rechits_FH+=(*rechit_energy)[i];
			    }
			  //.. Do something with each HGCAL rechit  ..//
			  h_hgcal_rechits->Fill((*rechit_energy)[i]);
			}
		      // else
		      // 	{ continue;}
		    } //noise cut
		} // channel 28 and 22 of chip -3
	    } // noisy flag
	} //Nrechits loop
      //cout<<"..."<<endl;
      //cout<<Esum_rechits_FH <<" "<<Esum_rechits_EE<<endl;
      h_EE_vs_FH->Fill(Esum_rechits_EE,Esum_rechits_FH);
      h_SEnergy_Rechits_EE->Fill(Esum_rechits_EE);
      h_SEnergy_Rechits_FH->Fill(Esum_rechits_FH);
      //    Esum_rechits_EE=0;
      Esum_rechits_FH=0;
      
      //////////////////////////////////////////
      
      
      //    cout<<NRechits<<"\t"<<count<<endl;
      
      ////////////////////////////////////////////
      //            AHCAL Part                  //
      ////////////////////////////////////////////
      
      /// Read AHCAL Tree
      if(DEBUG) cout<<"DEBUG: Start Analylizing AHCAL  RecHits!!"<<endl;
      // if(ahc_nHits>5)
      //   {
      // 	cout<<ahc_nHits<<" "<<runNumber<<" "<<jentry<<endl;}
      for(int i = 0 ; i < ahc_nHits; i++)
	{
	  //.. Do something with each AHCAL rechit  ..//
	  if((*ahc_hitEnergy)[i]>0.5)
	    {
	      
	      //	cout<<"..."<<(*rechit_layer)[i]<<endl;
	      //rechits_layer[(*rechit_layer)[i]+40]+=(*rechit_energy)[i];
	      Esum_rechits_AH+=(*ahc_hitEnergy)[i];
	      h_ahcal_rechits->Fill((*ahc_hitEnergy)[i]);
	    }
	}
      //cout<<Esum_rechits_AH<<endl;
      h_SEnergy_Rechits_AH->Fill(Esum_rechits_AH);
      h_EE_vs_AH->Fill(Esum_rechits_EE,Esum_rechits_AH);
      Esum_rechits_EE=0;
      Esum_rechits_AH=0;
      //loop over layers
      for(int ilayer=1;ilayer<41;ilayer++)
	{
	  //cout<<rechits_layer[ilayer]<<endl;
	  h_SEnergy_rechits_vs_layer->Fill(ilayer,rechits_layer[ilayer]);
	  rechits_layer[ilayer]=0;
	}
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      
      if(DEBUG && jentry > 50) return;
      
    
    } // loop over entries
  

  ///////////////////////////////////////////////////////////
  ///////  E N D     O F     E N T R Y     L O O P     //////
  ///////////////////////////////////////////////////////////
  
  cout<<"Got Out "<<jentry<<endl;


}

