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
      int Nrechits_FH_module_42=0,Nrechits_FH_module_45=0;
      for(int i =0;i<NRechits;i++)
	{
	  if(rechit_layer->at(i)==36)
	    {
	      if(rechit_module->at(i)==42 ) Nrechits_FH_module_42++;
	      else if(rechit_module->at(i)==45 ) Nrechits_FH_module_45++;
	      else continue;
	    }
	  else
	    continue;
	}
      // REJECTING EVENTS WITH VERY HIGH OCCUPANCY in FH9_P4 and FH9_P5 for config1      
      if(Nrechits_FH_module_42 >80 ||  Nrechits_FH_module_45 > 80) continue;

      Nrechits_FH_module_42=0;
      Nrechits_FH_module_45=0;


      ////////////////////////////////////////////
      //            HGCAL Part                  //
      ////////////////////////////////////////////
      
      /// Read HGCAL Tree
      if(DEBUG) cout<<"DEBUG: Start Analylizing HGCAL  RecHits!!"<<endl;
      // total_rechits_energy_EE=0;
      // total_rechits_energy_FH=0;
      ratio_inGeV_AH=0;Energy_EE_inGeV=0;
      Energy_FH_inGeV=0;
      Energy_AH_inGeV=0;
      totalEnergy_inGeV=0;
      ratio_inGeV_FH=0;
      ratio_AH_inMips=0;

      ////
      EnergySum_SSinEE=0;
      EnergySum_SSinEE_rmAHCAL=0;
      EnergySum_SSinFH=0;
      EnergySum_SSinFH_rmEE=0;
      EnergySum_MipLike=0;
      EnergySum_rmMipLike=0;
       EnergySum_SSinEE_inMips=0;
      EnergySum_SSinEE_rmAHCAL_inMips=0;
      EnergySum_SSinFH_inMips=0;
      EnergySum_SSinFH_rmEE_inMips=0;
      EnergySum_MipLike_inMips=0;
      EnergySum_rmMipLike_inMips=0;

      ///
      
      ratio_EE_inMips=0;
      ratio_FH_inMips=0;
      ratio_inGeV_EE=0;
      total_rechits_energy=0;
      total_rechits_energy_nocut=0;
      Esum_rechits_EE_allcut=0;
      Esum_rechits_FH_allcut=0;
      Esum_rechits_FH=0;
      Esum_rechits_EE=0;
      Esum_rechits_AH=0;
      Esum_rechits_AH_nocut=0;
      count_afterCuts=0;
      count=0;
      count_FH=0;
      count_FH_afterCuts=0;
      for( int i = 0 ; i < NRechits; i++)
      	{
      	  int channel = rechit_channel->at(i);
      	  int chip = rechit_chip->at(i);
      	  int en_chan = 1000*chip + channel;
      	  int layer = rechit_layer->at(i);
      	  float energy=rechit_energy->at(i);
      	  
	  if(rechit_noise_flag->at(i)==1) continue;      	       		
	  if(en_chan == 44 || en_chan == 3028 || en_chan == 3022) continue;		  		  
	  if(layer == 1 && chip == 0) continue;		  
	  if(rechit_module->at(i) == 39 && en_chan == 1016) continue;
	  if(rechit_energy->at(i) < 0.5) continue;
      	     
	  if(layer<=28)
            {
	      count_afterCuts++;
              Esum_rechits_EE_allcut+=energy;
              total_rechits_energy+=energy;
	    }
       
      	  else
      	    {
              count_FH_afterCuts++;
              Esum_rechits_FH_allcut+=energy;
      	      total_rechits_energy+=energy;
      	    }
      	  h_hgcal_rechits->Fill(rechit_energy->at(i));
      	} //Nrechits loop

      ratio= Esum_rechits_EE_allcut/Esum_rechits_FH_allcut;
      h_ratioEEvsFH->Fill(ratio);
      h_totalrechits_allcut->Fill(count_afterCuts);
      h_EE_vs_FH_allcut->Fill(Esum_rechits_EE_allcut,Esum_rechits_FH_allcut);
      ///energy in GeV
      Energy_EE_inGeV= 0.0105*Esum_rechits_EE_allcut;
      Energy_FH_inGeV= 0.0789*Esum_rechits_FH_allcut;
      h_Energy_EEvsFH_inGeV->Fill(Energy_EE_inGeV,Energy_FH_inGeV);
      h_Energy_EE_inGeV->Fill(Energy_EE_inGeV);

      h_Energy_FH_inGeV->Fill(Energy_FH_inGeV);

      h_SEnergy_Rechits_EE_allcut->Fill(Esum_rechits_EE_allcut);
      h_SEnergy_Rechits_FH_allcut->Fill(Esum_rechits_FH_allcut);
      
      
      // //////////////////////////////////////////
      
      
      // //    cout<<NRechits<<"\t"<<count<<endl;
      
      // ////////////////////////////////////////////
      // //            AHCAL Part                  //
      // ////////////////////////////////////////////
      
      /// Read AHCAL Tree
      if(DEBUG) cout<<"DEBUG: Start Analylizing AHCAL  RecHits!!"<<endl;

      for(int i = 0 ; i < ahc_nHits; i++) {
      	Esum_rechits_AH_nocut+=(*ahc_hitEnergy)[i];
      	//total_rechits_energy_nocut+=(*ahc_hitEnergy)[i];
      	if(ahc_hitK->at(i) == 38) continue;
	Esum_rechits_AH+=(*ahc_hitEnergy)[i];
	h_ahcal_rechits->Fill((*ahc_hitEnergy)[i]);
      	      //total_rechits_energy+=(*ahc_hitEnergy)[i];
      }//
      //cout<<Esum_rechits_AH<<endl;
      h_SEnergy_Rechits_AH->Fill(Esum_rechits_AH);
      h_EE_vs_AH->Fill(Esum_rechits_EE_allcut,Esum_rechits_AH);
      h_FH_vs_AH->Fill(Esum_rechits_FH_allcut,Esum_rechits_AH);
      // h_EE_vs_AH_nocut->Fill(Esum_rechits_EE_nocut,Esum_rechits_AH_nocut);
      // h_FH_vs_AH_nocut->Fill(Esum_rechits_FH_nocut,Esum_rechits_AH_nocut);
      ///in GeV
      Energy_AH_inGeV=0.0316*Esum_rechits_AH;
      h_Energy_AH_inGeV->Fill(Energy_AH_inGeV);
      
      h_Energy_EEvsAH_inGeV->Fill(Energy_EE_inGeV,Energy_AH_inGeV);
      h_Energy_FHvsAH_inGeV->Fill(Energy_FH_inGeV,Energy_AH_inGeV);
      h_total_rechitsEenrgy->Fill(total_rechits_energy);
      totalEnergy_inGeV= Energy_EE_inGeV + Energy_FH_inGeV; //+ Energy_AH_inGeV;
      ratio_inGeV_AH=Energy_AH_inGeV/totalEnergy_inGeV;
      ratio_AH_inMips=Esum_rechits_AH/total_rechits_energy;
      h_ratioInGeV_AH->Fill(ratio_inGeV_AH);
      h_ratioInMips_AH->Fill(ratio_AH_inMips);
      h_totalEnergy_inGeV->Fill(totalEnergy_inGeV);
      ratio_inGeV_FH = Energy_FH_inGeV/(totalEnergy_inGeV);
      ratio_inGeV_EE = Energy_EE_inGeV/(totalEnergy_inGeV);
      h_ratioInGeV_EE->Fill(ratio_inGeV_EE);
      h_ratioInGeV_FH->Fill(ratio_inGeV_FH);
      //check if it is a mip like or not ///
      /////           condition for MIP-LIke   ///
      ////  if mips in EE <100 and in FH <60

      if(Esum_rechits_EE_allcut>100)// &&  Esum_rechits_FH_allcut <60) /// SS in EE
	{
	  //cout<<Esum_rechits_FH_allcut<<endl;
	  EnergySum_SSinEE = 0.0105*Esum_rechits_EE_allcut + 0.0789*Esum_rechits_FH_allcut+0.0316*Esum_rechits_AH;
	  EnergySum_SSinEE_rmAHCAL = 0.0105*Esum_rechits_EE_allcut + 0.0789*Esum_rechits_FH_allcut;
	  EnergySum_rmMipLike= 0.0105*Esum_rechits_EE_allcut + 0.0789*Esum_rechits_FH_allcut+0.0316*Esum_rechits_AH;

	  EnergySum_SSinEE_inMips = Esum_rechits_EE_allcut + Esum_rechits_FH_allcut+Esum_rechits_AH;
          EnergySum_SSinEE_rmAHCAL_inMips= Esum_rechits_EE_allcut + Esum_rechits_FH_allcut;
          EnergySum_rmMipLike_inMips= Esum_rechits_EE_allcut + Esum_rechits_FH_allcut+Esum_rechits_AH;
	  ///
	  h_EnergySum_SSinEE->Fill(EnergySum_SSinEE);
	  h_EnergySum_SSinEE_rmAHCAL->Fill(EnergySum_SSinEE_rmAHCAL);
	  h_EnergySum_rmMipLike->Fill(EnergySum_rmMipLike);

	  h_EnergySum_SSinEE_inMips->Fill(EnergySum_SSinEE_inMips);
	  h_EnergySum_SSinEE_rmAHCAL_inMips->Fill(EnergySum_SSinEE_rmAHCAL_inMips);
	  h_EnergySum_rmMipLike_inMips->Fill(EnergySum_rmMipLike_inMips);


	}
      else if( Esum_rechits_EE_allcut<100 && Esum_rechits_FH_allcut >60) // SS in FH
	{
	  EnergySum_SSinFH = 0.0105*Esum_rechits_EE_allcut + 0.0789*Esum_rechits_FH_allcut+0.0316*Esum_rechits_AH;
	  EnergySum_SSinFH_rmEE =  0.0789*Esum_rechits_FH_allcut+0.0316*Esum_rechits_AH;
	  EnergySum_rmMipLike= 0.0105*Esum_rechits_EE_allcut + 0.0789*Esum_rechits_FH_allcut+0.0316*Esum_rechits_AH;

	  EnergySum_SSinFH_inMips = Esum_rechits_EE_allcut + Esum_rechits_FH_allcut+Esum_rechits_AH;
          EnergySum_SSinFH_rmEE_inMips=  Esum_rechits_FH_allcut+Esum_rechits_AH;
          EnergySum_rmMipLike_inMips= Esum_rechits_EE_allcut + Esum_rechits_FH_allcut+Esum_rechits_AH;

	  h_EnergySum_SSinFH->Fill( EnergySum_SSinFH);
	  h_EnergySum_SSinFH_rmEE->Fill( EnergySum_SSinFH_rmEE);
	  h_EnergySum_rmMipLike->Fill(EnergySum_rmMipLike);
	  h_EnergySum_rmMipLike_inMips->Fill(EnergySum_rmMipLike_inMips);
	  h_EnergySum_SSinFH_inMips->Fill( EnergySum_SSinFH_inMips);
	  h_EnergySum_SSinFH_rmEE_inMips->Fill( EnergySum_SSinFH_rmEE_inMips);

	}
      else // mip like
	{
	  EnergySum_MipLike = 0.0105*Esum_rechits_EE_allcut + 0.0789*Esum_rechits_FH_allcut+0.0316*Esum_rechits_AH;
	  EnergySum_MipLike_inMips = Esum_rechits_EE_allcut + Esum_rechits_FH_allcut+Esum_rechits_AH;
	  h_EnergySum_MipLike->Fill(EnergySum_MipLike);
	  h_EnergySum_MipLike_inMips->Fill(EnergySum_MipLike_inMips);
	}

      // h_EnergySum_SSinEE->Fill(EnergySum_SSinEE);
      // h_EnergySum_SSinEE_rmAHCAL->Fill(EnergySum_SSinEE_rmAHCAL);
      // h_EnergySum_SSinFH->Fill( EnergySum_SSinFH);
      // h_EnergySum_SSinFH_rmEE->Fill( EnergySum_SSinFH_rmEE);
      // h_EnergySum_MipLike->Fill(EnergySum_MipLike);
      // h_EnergySum_rmMipLike->Fill(EnergySum_rmMipLike);

      //  h_EnergySum_SSinEE_inMips->Fill(EnergySum_SSinEE_inMips);
      // h_EnergySum_SSinEE_rmAHCAL_inMips->Fill(EnergySum_SSinEE_rmAHCAL_inMips);
      // h_EnergySum_SSinFH_inMips->Fill( EnergySum_SSinFH_inMips);
      // h_EnergySum_SSinFH_rmEE_inMips->Fill( EnergySum_SSinFH_rmEE_inMips);
      // h_EnergySum_MipLike_inMips->Fill(EnergySum_MipLike_inMips);
      // h_EnergySum_rmMipLike_inMips->Fill(EnergySum_rmMipLike_inMips);

      /// an additional cut          ///
      // for positrons, events that has >100 mIps in AHCAL, reject those//

      
      // if(Esum_rechits_AH>1000)
      // 	{
      // 	  cout<<runNumber<<"\t"<<eventNumber<<endl;
      // 	}
      //in mips                                                                                                                                                           
      ratio_EE_inMips=Esum_rechits_EE_allcut/total_rechits_energy;
      ratio_FH_inMips=Esum_rechits_FH_allcut/total_rechits_energy;
      h_ratioInMips_EE->Fill(ratio_EE_inMips);
      h_ratioInMips_FH->Fill(ratio_FH_inMips);

            // h_total_rechitsEenrgy_nocut->Fill(total_rechits_energy_nocut);
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      
      if(DEBUG && jentry > 50) return;
      //      gSystem->Exit(0);
      //if(jentry > 10000) break;
    
    } // loop over entries
  

  ///////////////////////////////////////////////////////////
  ///////  E N D     O F     E N T R Y     L O O P     //////
  ///////////////////////////////////////////////////////////
  //  gSystem->Exit(0);
  cout<<"Got Out "<<jentry<<endl;


}

