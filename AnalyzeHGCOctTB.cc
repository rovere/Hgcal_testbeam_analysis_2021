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

  if (argc < 3) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "energy" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *energy = argv[3];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, energy);
  //cout << "dataset " << data << " " << endl;
  //cout << "configuration " << config << " " << endl;
  cout << "energy " << energy << " " << endl;
  // int chi2_method = atoi(energy);
  hgcOctTB.EventLoop(energy);//, min_, max_);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *energy) {


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries3 = fChain3->GetEntriesFast();
  Long64_t hgc_jentry = 0;
  //  int Elist[85] =  {200, 14 , 18 , 22 , 26 , 30,  34 , 38 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
  // 		    82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
  // 		    154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
  // 		    226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
  // 		    298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};
  float Elist1[9]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};//,120,200,250,300};
  int Elist[8]={20,50,80,100,120,200,250,300};
  cout << "nentries " << nentries << endl;
  //  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  Long64_t nbytes3 = 0, nb3 = 0;
  int count_evt[9]={};
  int decade = 0;
  
  bool DEBUG=false;
  //counter
  int nHgc=0, nAhc=0, nrechits=0;
  Long64_t jentry;
  int chi2_method = atoi(energy);
  int count[350]={};
  float lambda[79];
  // nentries=8000;
  float ee_rescaling = 1.035;
  float fh_rescaling = 1.095;
  float ah_rescaling = 1.095;

  // chi2 method
  std::vector<std::vector<ROOT::Math::SVector<double, 3> > > consts_EH_vec_bin;
  std::vector<std::vector<ROOT::Math::SVector<double, 3> > > values_EH_vec_bin;
  std::vector<std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > >coeff_EH_vec_bin;
                                                                          
  std::vector<std::vector<ROOT::Math::SVector<double, 2> > > consts_H_vec_bin;
  std::vector<std::vector<ROOT::Math::SVector<double, 2> > > values_H_vec_bin;
  std::vector<std::vector<ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2>> > >coeff_H_vec_bin;
  bool isInverted_EH = false;
  bool isInverted_H = false;

  for (int j =0 ;j<8;j++)
    {
  std::vector<ROOT::Math::SVector<double, 3> >consts_EH_vec;
  std::vector<ROOT::Math::SVector<double, 3> >values_EH_vec;
  std::vector<ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3>> > coeff_EH_vec;
  ROOT::Math::SMatrix<double,3, 3, ROOT::Math::MatRepStd<double,3> > coeffs_EH;
  ROOT::Math::SVector<double, 3> consts_EH;
  ROOT::Math::SVector<double, 3> values_EH;

  //  bool isInverted_EH = false;
  for(int i =0; i<9; i++)
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
  consts_EH_vec_bin.push_back(consts_EH_vec);
  values_EH_vec_bin.push_back(values_EH_vec);
  coeff_EH_vec_bin.push_back(coeff_EH_vec);

  std::vector<ROOT::Math::SVector<double, 2> >consts_H_vec;
  std::vector<ROOT::Math::SVector<double, 2> >values_H_vec;
  std::vector<ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2>> > coeff_H_vec;
  bool isInverted_H = false;
  ROOT::Math::SMatrix<double,2, 2, ROOT::Math::MatRepStd<double,2> > coeffs_H;
  ROOT::Math::SVector<double, 2> consts_H;
  ROOT::Math::SVector<double, 2> values_H;

  for(int i =0; i<9; i++)
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
  consts_H_vec_bin.push_back(consts_H_vec);
  values_H_vec_bin.push_back(values_H_vec);
  coeff_H_vec_bin.push_back(coeff_H_vec);

    }
  int ahcal_layer[10]= {43,47,51,55,59,63,67,71,75,79};
  for (jentry=0; jentry<nentries;jentry++,hgc_jentry++)
  {
    //   for (jentry=0; jentry<10000;jentry++,hgc_jentry++) {
    // ==============print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;
      // if(jentry%10000==0)
	
      // 	cout<<"entre:  "<<jentry<<endl;
      // ===============read this entry == == == == == == == == == == ==

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { break; cout<<"Breaking"<<endl;}
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //cout<<"insides the event loop: check1"<<endl;

      trueBeamEnergy = beamEnergy;
      //if(beamEnergy!=200) continue;
      ////   MESSAGE ////
      // apparently code is not running beyond this point ///
      int nrechits[5]={};
      int Nrechit_en[5][8]={};
      float total_Energy[5][8]={};
      double dr_cutValues[5]={0.5,0.1,0.2,0.3,0.5};
      int total_nrechits[5][8]={};
      float total_Energy_EE[5][8]={};
      float total_Energy_FH[5][8]={};
      float total_Energy_AH[5][8]={};
      float total_Energy_EE_mips[5][8]={};
      float total_Energy_FH_mips[5][8]={};
      float total_Energy_AH_mips[5][8]={};
      float total_Energy_mips[5][8]={};
      float total_EE_mips[5][8]={};
      float total_FH_mips[5][8]={};
      float total_AH_mips[5][8]={};
      float EnergyRechits_EE =0.0,EnergyRechits_FH=0.0,EnergyRechits_AH=0.0;
      //      cout<<"before HGCAL hits loop"<<endl;
      for(int i = 0 ; i < rechit_channel->size(); i++){
	int channel = rechit_channel->at(i);
	int chip = rechit_chip->at(i);
	int en_chan = 1000*chip + channel;
	int layer = rechit_layer->at(i);
	float energy=rechit_energy->at(i);
	if(layer<=28)

	  EnergyRechits_EE+=(energy/ee_rescaling);
	else
	  EnergyRechits_FH+=(energy/fh_rescaling);
	// if(layer<=28)
	//   total_EE_mips+=energy;
	// else
	//   total_FH_mips+=energy;
	for(int j=0;j<1;j++)
	  {
	    nrechits[j]++;
	    h_rechitEnergy_Mips[j]->Fill(energy);
	    h_rechitX[j]->Fill(rechit_x->at(i));
	    h_rechitY[j]->Fill(rechit_y->at(i));
	    h_rechitZ[j]->Fill(rechit_z->at(i));
	    h_rechitXvsY[j]->Fill(rechit_x->at(i),rechit_y->at(i));
	    h_rechitXvsZ[j]->Fill(rechit_z->at(i),rechit_x->at(i));		  
	    for (int k=0;k<8;k++)
	      {
		if((Elist[k]-2)<=trueBeamEnergy && trueBeamEnergy<=(Elist[k]+2))
		  {

		    if(layer<=28)
		      {
		      total_EE_mips[j][k]+=(energy/ee_rescaling);
		      total_Energy_mips[j][k]+=(energy/ee_rescaling);
		      }
		    else
		      {
		      total_FH_mips[j][k]+=(energy/fh_rescaling);
		      total_Energy_mips[j][k]+=(energy/fh_rescaling);
		      }

		    //  cout<<"genEt"<<genEt<<"\t"<<k<<"\t"<<Elist[k]<<endl;
		    total_nrechits[j][k]++;
		    //		    total_Energy_mips[j][k]+=energy;
		    h_bin_rechitEn_Mips[j][k]->Fill(energy);
		    h_bin_rechitX[j][k]->Fill(rechit_x->at(i));
		    h_bin_rechitY[j][k]->Fill(rechit_y->at(i));
		    h_bin_rechitZ[j][k]->Fill(rechit_z->at(i));

		    //		    total_Energy[j][k]+=(0.001*energy[i]);
		    if(layer<=28)
		      {
			h_bin_rechitEn_Mips_EE[j][k]->Fill(energy);
			//h_bin_rechitEn_EE[j][k]->Fill(0.001*hit_en[i]);
			h_bin_rechitX_EE[j][k]->Fill(rechit_x->at(i));
			h_bin_rechitY_EE[j][k]->Fill(rechit_y->at(i));
			h_bin_rechitZ_EE[j][k]->Fill(rechit_z->at(i));
			//	total_Energy_EE[j][k]+=(0.001*energy[i]);
			//total_Energy_EE_mips[j][k]+=energy;
			h_bin_rechitXvsY_EE[j][k]->Fill(rechit_x->at(i),rechit_y->at(i));
			h_bin_rechitXvsZ_EE[j][k]->Fill(rechit_z->at(i),rechit_x->at(i));
			//total_EE_mips+=energy;
		      }
		    else //if(hit_lay[i]<=38 && hit_lay[i]>26)
		      {
			//total_FH_mips+=energy;
			h_bin_rechitEn_Mips_FH[j][k]->Fill(energy);
			//h_bin_rechitEn_FH[j][k]->Fill(0.001*hit_en[i]);
			h_bin_rechitX_FH[j][k]->Fill(rechit_x->at(i));
			h_bin_rechitY_FH[j][k]->Fill(rechit_y->at(i));
			h_bin_rechitZ_FH[j][k]->Fill(rechit_z->at(i));
			//total_Energy_FH[j][k]+=(0.001*hit_en[i]);
			//total_Energy_FH_mips[j][k]+=energy;
			h_bin_rechitXvsY_FH[j][k]->Fill(rechit_x->at(i),rechit_y->at(i));
			h_bin_rechitXvsZ_FH[j][k]->Fill(rechit_z->at(i),rechit_x->at(i));

		      }
		  }
	      }
	  }
      }

      //cout<<"after HGCAL hits loop"<<endl;
      for(int i = 0 ; i <ahc_hitEnergy->size(); i++) {
	int layer = ahc_hitK->at(i)+40;
	for(int i_l=0; i_l<10; i_l++)
	  {
	    if(layer ==ahcal_layer[i_l])
	      {
		//total_AH_mips+=ahc_hitEnergy->at(i);
		EnergyRechits_AH+=(ahc_hitEnergy->at(i)/ah_rescaling);
		for(int j=0;j<1;j++)
		  {
		    h_rechitEnergy_Mips[j]->Fill(ahc_hitEnergy->at(i));
		    h_rechitX[j]->Fill(ahc_hitX->at(i));
		    h_rechitY[j]->Fill(ahc_hitY->at(i));
		    h_rechitZ[j]->Fill(ahc_hitZ->at(i));
		    h_rechitXvsY[j]->Fill(ahc_hitZ->at(i),ahc_hitX->at(i));
		    //total_Energy_AH_mips+=ahc_hitEnergy->at(i);
		    for (int k=0;k<8;k++)
		      {
			if((Elist[k]-2)<=trueBeamEnergy && trueBeamEnergy<=(Elist[k]+2))
			  {
			    total_AH_mips[j][k]+=(ahc_hitEnergy->at(i)/ah_rescaling);
			    total_nrechits[j][k]++;
			    total_Energy_mips[j][k]+=(ahc_hitEnergy->at(i)/ah_rescaling);
			    h_bin_rechitEn_Mips[j][k]->Fill(ahc_hitEnergy->at(i));
			    h_bin_rechitX[j][k]->Fill(ahc_hitX->at(i));
			    h_bin_rechitY[j][k]->Fill(ahc_hitY->at(i));
			    h_bin_rechitZ[j][k]->Fill(ahc_hitZ->at(i));
			    
		    //total_Energy_mips[j][k]+=ahc_hitEnergy->at(i);
		    
		    h_bin_rechitEn_Mips_AH[j][k]->Fill(ahc_hitEnergy->at(i));
		    //h_bin_rechitEn_AH[j][k]->Fill(0.001*hit_en[i]);//
		    h_bin_rechitX_AH[j][k]->Fill(ahc_hitX->at(i));
		    h_bin_rechitY_AH[j][k]->Fill(ahc_hitY->at(i));
		    h_bin_rechitZ_AH[j][k]->Fill(ahc_hitZ->at(i));
		    //total_Energy_AH_mips[j][k]+=ahc_hitEnergy->at(i);//(0.001*hit_en[i]);
		    //		    total_Ahc_HitEnergy->At(I)_AH_mips[j][k]+=ahc_hitEnergy->at(i);
		    h_bin_rechitXvsY_AH[j][k]->Fill(ahc_hitX->at(i),ahc_hitY->at(i));
		    h_bin_rechitXvsZ_AH[j][k]->Fill(ahc_hitZ->at(i),ahc_hitX->at(i));
			      
		  }
		
	      }
		    
	  }
	      }		  
	  }
      }
      
      float EnergyEE_GeV=0.0,EnergyFH_GeV=0.0,EnergyAH_GeV=0.0,total_E_GeV=0.0;
      EnergyEE_GeV = 0.0105*EnergyRechits_EE;
      EnergyFH_GeV = 0.0812*EnergyRechits_FH;
      EnergyAH_GeV = 0.12504*EnergyRechits_AH;
      double EE_detscale = EnergyEE_GeV; 
      double FH_detscale = EnergyFH_GeV;//(rechitEnergySum_FH/FH_AH_scale);                                                              
      double AH_detscale = EnergyAH_GeV ;//Esum_rechits_AH_inGeV;//(alpha_*rechitEnergySum_AH)/FH_AH_scale; 
      double full_energy = FH_detscale+AH_detscale;
      double total_energy = EE_detscale + FH_detscale + AH_detscale;
      total_E_GeV = EE_detscale + FH_detscale + AH_detscale;
      double E_beam =0.0;
      double O = 0.0;
      double sigma2_H = (0.09*0.09)*(total_energy*total_energy) + (1.228*1.228)*total_energy;
      double tot_E_gev= total_energy;
      for (int j =0; j<8;j++)
	{
	  if(trueBeamEnergy==Elist[j])
	    {
	      h_rechitEE_vsHad_GeV[j]->Fill(EnergyEE_GeV, EnergyFH_GeV+EnergyAH_GeV);
	      hprof_EEvsHad[j]->Fill(EnergyEE_GeV,EnergyFH_GeV+EnergyAH_GeV);
	      if(rechit_shower_start_layer>28)
		h_rechitEE_vsHad_GeV_MipsInEE[j]->Fill(EnergyEE_GeV,EnergyFH_GeV+EnergyAH_GeV);
	      for(int i_bin=0;i_bin<9;i_bin++)
		{
		  int xmax = 0;
		  if(i_bin==0)
		    xmax =0;
		  else
		    xmax = Elist1[i_bin-1]*Elist[j];
		  
		  if(xmax<EE_detscale && EE_detscale<=(Elist[j]*Elist1[i_bin]))
		{
		  //		  count_evt[i_bin]+=1;
		  //cout<<EE_detscale<<endl;
		  E_beam = Elist[j];
		  
		  // if(rechit_shower_start_layer<=28)
		  //   {
		  //     float w1 = getChi2Weights_EH(Elist1[i_bin]+2).at(0);
		  //     float w2 = getChi2Weights_EH(Elist1[i_bin]+2).at(1);
		  //     float w3 = getChi2Weights_EH(Elist1[i_bin]+2).at(2);
		  //     //cout<<w1<<"\t"<<w2<<"\t"<<w3<<endl;
		  //     // h_totalSum_EE[i]->Fill(w1*EE_detscale);
		  //     // h_totalSum_FH[i]->Fill(w2*FH_detscale);
		  //     // h_totalSum_AH[i]->Fill(w3*AH_detscale);
		  //     h_totalSum_total[i]->Fill(w1*EE_detscale+w2*FH_detscale+w3*AH_detscale);
		  //   }
		  if(rechit_shower_start_layer>28)
		    {
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
      	      coeffs_H(0,0) = (FH_detscale*FH_detscale)/sigma2_H;
      	      coeffs_H(0,1) = (AH_detscale*FH_detscale)/sigma2_H;
      	      coeffs_H(1,0) = (FH_detscale*AH_detscale)/sigma2_H;
      	      coeffs_H(1,1) = (AH_detscale*AH_detscale)/sigma2_H;

      	      consts_H(0)   = (E_beam-O)*FH_detscale/sigma2_H;
      	      consts_H(1)   = (E_beam-O)*AH_detscale/sigma2_H;
      	      consts_H_vec_bin[j][i_bin] += consts_H;
      	      coeff_H_vec_bin[j][i_bin]+=coeffs_H;
      	    }
      	  else
      	    {
      	      double sigma2_EH = (0.084*0.084) + (1.388*1.388)*tot_E_gev;
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
      	      consts_EH_vec_bin[j][i_bin]+= consts_EH;
      	      coeff_EH_vec_bin[j][i_bin]+=coeffs_EH;
      	    }
		}
      	}
	}
  
	}
      // //cout<<"After ACAL hits loop"<<endl;	
      // // 	    }
      // // 	}
      // // //cout<<nhits<<endl;
      // // //h_nhits->Fill(nhits);
      // //      h_bin_rechitEn_EEvsFH->Fill(total_EE_mips[5],(total_FH_mips[5]+1.54*total_AH_mips[5]));
      // //      cout<<
      // //if(198<=trueBeamEnergy&& trueBeamEnergy<=202)
      // //	h_bin_rechitEn_EEvsFH->Fill(total_EE_mips[5],(total_FH_mips[5]+(1.54*total_AH_mips[5])));
      for(int j=0;j<1;j++)
      	{
      	  h_nhits[j]->Fill(nrechits[j]);
      	  for (int k=0;k<8;k++)
      	    {
      	      if(Elist[k]-2<=trueBeamEnergy&& trueBeamEnergy<=Elist[k]+2)
      	  	{
		  //if(k==5)
		    //		    h_bin_rechitEn_EEvsFH->Fill(total_EE_mips[5],(total_FH_mips[5]+1.54*total_AH_mips[5])); 
		  
		  //h_bin_rechitEn_EEvsFH->Fill(total_EE_mips[5],(total_FH_mips[5]+(1.54*total_AH_mips[5])));
		  h_bin_rechitEn_EEvsFH[j][k]->Fill(total_EE_mips[j][k],(total_FH_mips[j][k]+1.54*total_AH_mips[j][k]));
      	  	  h_bin_TotalrechitEn[j][k]->Fill( total_Energy[j][k]);
      		  h_bin_TotalrechitEn_EE[j][k]->Fill( total_Energy_EE[j][k]);
      		  h_bin_TotalrechitEn_FH[j][k]->Fill( total_Energy_FH[j][k]);
      		  h_bin_TotalrechitEn_AH[j][k]->Fill( total_Energy_AH[j][k]);
		  
      		  float ratio= total_Energy[j][k]/trueBeamEnergy;
      		  float ratio_EE= total_Energy_EE[j][k]/trueBeamEnergy;
      		  float ratio_FH= total_Energy_FH[j][k]/trueBeamEnergy;
      		  float ratio_AH= total_Energy_AH[j][k]/trueBeamEnergy;
      		  h_bin_ratiorechitEn[j][k]->Fill(ratio);
                  h_bin_ratiorechitEn_EE[j][k]->Fill(ratio_EE);
                  h_bin_ratiorechitEn_FH[j][k]->Fill(ratio_FH);
                  h_bin_ratiorechitEn_AH[j][k]->Fill(ratio_AH);

		  h_bin_total_Mips[j][k]->Fill( total_Energy_mips[j][k]);
		  h_bin_total_Mips_EE[j][k]->Fill( total_Energy_EE_mips[j][k]);
		  h_bin_total_Mips_FH[j][k]->Fill( total_Energy_FH_mips[j][k]);
		  h_bin_total_Mips_AH[j][k]->Fill( total_Energy_AH_mips[j][k]);

      	  	  h_bin_nrechits[j][k]->Fill(total_nrechits[j][k]);
		  
      	  	}
      	    }
      	}
      total_Energy[5][8]={};
      total_Energy_EE[5][8]={};
      total_Energy_FH[5][8]={};
      total_Energy_AH[5][8]={};
      total_Energy_mips[5][8]={};
      total_Energy_EE_mips[5][8]={};
      total_Energy_FH_mips[5][8]={};
      total_Energy_AH_mips[5][8]={};


      // totalEnergy_inGeV=0;
      // EnergySum_SSinEE=0;
      // EnergySum_SSinFH=0;
      // Esum_rechits_FH=0;
      // Esum_rechits_EE=0;
      // Esum_rechits_AH=0;
      // Esum_rechits_FH_inGeV=0;
      // Esum_rechits_EE_inGeV=0;
      // Esum_rechits_AH_inGeV=0;
      
      // ////////////////////////////////////////////
      // //            HGCAL Part                  //
      // ////////////////////////////////////////////
      // total_rechits=0;
      // rechits_EE=0;
      // rechits_FH=0;
      // rechits_AH=0;
      // if(event ==3482)
      // 	{	//	continue;
      if(DEBUG) cout<<"DEBUG: Start Analylizing HGCAL  RecHits!!"<<endl;
      // total_rechits_energy_EE=0;
      // total_rechits_energy_FH=0;
      if(DEBUG) cout<<"DEBUG: End of Entry = "<<jentry<<endl;
      if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
      
      //if(DEBUG && jentry > 10) break;
      //      gSystem->Exit(0);
      //if(jentry > 2) break;
      //pion_tree->Fill();
    } // loop overentries
  //  char* name = new char[1000];


  // cout<<endl<<endl;


  // ///////////////////////////////////////////////////////////
  // ///////  E N D     O F     E N T R Y     L O O P     //////
  // ///////////////////////////////////////////////////////////
  // //  gSystem->Exit(0);
  char* name1 = new char[1000];
  sprintf(name1,"./events_vs_energy.txt");
  std::ofstream file_H_;
  file_H_.open(name1,ios::out);

  cout<<"Got Out "<<endl;
  for (int i=10;i<351;i++)
    {
      file_H_<<i<<"\t"<<count[i]<<"\n";
    }

  for (int i_j=0;i_j<8;i_j++)
    {
      char* name = new char[1000];
      sprintf(name,"chi2_calibFact_H_hadrons_flatEn_%d_%d_GeV.txt",chi2_method,Elist[i_j]);
      std::ofstream file_H;
      file_H.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_H_hadrons_flatEn_%d.txt",chi2_method,ios::out);                               
      sprintf(name,"chi2_calibFact_EH_hadrons_flatEn_%d_%d_GeV.txt",chi2_method,Elist[i_j]);
      std::ofstream file_EH;
      file_EH.open(name,ios::out);//"../txt_maps/chi2_map/chi2_calibFact_EH_hadrons_flatEn_%d.txt",chi2_method,ios::out  
      
      for(int i_en =0; i_en<9; i_en++)
	{
      cout<<"count :   "<<count_evt[i_en]<<endl;
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
      consts_EH = consts_EH_vec_bin[i_j][i_en];
      values_EH = values_EH_vec_bin[i_j][i_en];
      coeffs_EH = coeff_EH_vec_bin[i_j][i_en];
      isInverted_EH = coeffs_EH.Invert();
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          cout<<coeffs_EH(i, j)<<"\t";
        }
        cout<<endl;
      }

      cout<<" "<<consts_EH(0)<<"\t"<<consts_EH(1)<<"\t"<<consts_EH(2)<<endl;
      consts_H = consts_H_vec_bin[i_j][i_en];
      values_H = values_H_vec_bin[i_j][i_en];
      coeffs_H = coeff_H_vec_bin[i_j][i_en];

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
        cout<<"EH Hadrons => For E = "<<Elist[i_en]<<"GeV, w1 = "<<values_EH(0)<<" ;w2 = "<<values_EH(1)<<" ;w3 = "<<values_EH(2)<<";"<<endl;
        file_EH<<Elist1[i_en]<<"\t"<<values_EH(0)<<"\t"<<values_EH(1)<<"\t"<<values_EH(2)<<"\n";
      }
      else {
        cout<<"Error: Could not invert for EH hadrons..."<<endl;
      }
      //isInverted_H = true;                                                                                                                      
      if(isInverted_H) {
        values_H = coeffs_H*consts_H;
        //cout<<"H Hadrons => For E = "<<E_beam<<"GeV, w1 = "<<values_H(0)<<" ;w2 = "<<values_H(1)<<" ;w3 = "<<values_H(2)<<";"<<endl;            
        cout<<"H Hadrons => For E = "<<Elist[i_en]<<"GeV, w1 = 0.0 ;w2 = "<<values_H(0)<<" ;w3 = "<<values_H(1)<<";"<<endl;
        file_H<<Elist1[i_en]<<"\t"<<values_H(0)<<"\t"<<values_H(1)<<"\n";
      }
      else {
        cout<<"Error: Could not invert for H hadrons..."<<endl;
      }
	}
    }
  // cout<<count_fh<<endl;
  // cout<<count_badtrack<<endl;
  // oFile->cd();
  // oFile->Write();
  // oFile->Close();

}
