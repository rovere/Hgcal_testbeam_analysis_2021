#ifndef AnalyzeHGCOctTB_H
#define AnalyzeHGCOctTB_H
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

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
  AnalyzeHGCOctTB(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *config="alpha", const char* chi_method = "-1", const char* folderID = "-1");
  ~AnalyzeHGCOctTB();
  //Bool_t   FillChain(TChain *chain, TChain *chain2, TChain *chain3, const TString &inputFileList);
   Bool_t   FillChain(TChain *chain, const TString &inputFileList); 
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *, const char *, const char *);
  void     BookHistogram(const char *, const char *, const char* chi_method);

  void moduleMap_init(const char *);
  void Alignment_Map_Init();
  void Noise_Map_Init();
  void layerPosition_init();
  void Chi2_Weight_Map_Init(int chi2_method);
  
  std::vector<bool> *noise_flag;
  TFile *oFile;
  const char *conf_;  
  int inEnergy_;
  int event_count[7]={};
  int count=0, count_afterCuts =0;//,count_nflag=0,count_nchannel=0,count_nchanLayer1=0, count_nflag0=0;
  float  totalEnergy_inGeV=0;
  float  EnergySum_SSinEE=0;
  float  EnergySum_SSinFH=0;
  float Esum_rechits_FH=0;
  float Esum_rechits_EE=0;
  float Esum_rechits_AH=0;
  float Esum_rechits_FH_inGeV=0;
  float Esum_rechits_EE_inGeV=0;
  float Esum_rechits_AH_inGeV=0;
  int total_rechits=0, rechits_EE=0, rechits_FH=0, rechits_AH=0;
  
  int count_FH=0,count_fh=0, count_badtrack=0;//,count_adiMask=0;
  int count_FH_afterCuts=0;
  TH1F* h_beamenergy;
  TH1F* h_particle;

  TH1F* h_true_beamenergy[6];
  TH1F* h_EnergySum_inEE;
  TH1F* h_EnergySum_inFH;
  TH1F* h_EnergySum_inAH;
  TH1F* h_EnergySum_inAH_inGeV;
  TH1F* h_EnergySum_inFH_inGeV;
  TH1F* h_EnergySum_inEE_inGeV;
  TH1F* h_EnergySum_inGeV;
  TH1F* h_EnergySum_ratio_inGeV;
  TH1F* h_EnergySum_ratio_inMips;
  TH1F* h_EnergySum_Rechit;
  TH1F* h_EnergySum_Rechit_normalized;
  TH2F* h_EnergySum_inGeVvstrue;
  TH2F* h_EnergySum_inMipsvstrue;
  TH2F* h_EnergySum_inAH_vs_true;
  TH2F* h_EnergySum_inFH_vs_true;
  TH2F* h_EnergySum_inEE_vs_true;
  TH2F* h_EnergySum_inAH_Mipsvs_true;// h_EnergySum_inAHMipsvs_true;
  TH2F* h_EnergySum_inFH_Mipsvs_true;
  TH2F* h_EnergySum_inEE_Mipsvs_true;
  TH1F* h_nrechits;
  TH1F* h_nrechits_EE;
  TH1F* h_nrechits_FH;
  TH1F* h_nrechits_AH;
  TH1F* h_ratio_nrechits_AH;
  TH1F* h_ratio_nrechits_FH;
  TH1F* h_ratio_nrechits_EE;
  
   TH1F* h_rechits_mips_EE;
  TH1F* h_rechits_mips_FH;
  TH1F* h_rechits_mips_AH;
  TH2F* h_rechitvslambda;  
  TH1F* h_rechits_GeV_EE;
  TH1F* h_rechits_GeV_FH;
  TH1F* h_rechits_GeV_AH;

  TH2F* h_nrechits_vs_EE;
  TH2F* h_nrechits_vs_FH;
  TH2F* h_nrechits_vs_AH;
  TH2F* h_EnergySum_inEE_vs_FH;
  TH2F* h_EnergySum_inEE_vs_FH_inGeV;
  TH2F* h_nrechits_EE_vs_FH;
  TH2F* h_EnergySum_ratio_vs_true;
  TH1F* h_EnergySum_ratio_flipped;
  TH2F*  h_EnergySum_ratio_vs_true_flipped;
  TH1F* h_ssinEE;
  TH1F* h_mipsinEE;
  TDirectory* resp;
  TDirectory* fix_wten;
  TH1F* hist_fixwt_total[85];
  TH1F* hist_fixwt_SS_EE[85];
  TH1F* hist_fixwt_SS_FH[85];

  TDirectory* chi2_mehtod1;
  TH1F* h_rechit_energy_FB_rel_weightScan[50];
  TH2F* h_muo_dxvsdy;
  TH1F* h_dR;
  TH2F* h_muo_dxvsdy_afterCorr;
  TH1F* h_dR_afterCorr;

  // 2d histograms for each energy bin
  //histograms for muons
  TH2F* h_xVs_trackX[3]; //at CE-E 1st, CE-H 2nd, AH 3rd
  TH2F* h_xVs_1trackX[3];
  TH2F* h_yVs_trackY[3];
  TH2F* h_yVs_1trackY[3]; //done
  //h_xVs_trackX  
  TH2F* h_dXvsdY_EE[2]; //EE = 1,14
  TH2F* h_dXvsdY_FH[2];  //FH = 1,6
  TH2F* h_dXvsdY_AH[2]; //AH 1, 20
  TH2F* h_dXvsdY_EE_Acorr[2]; //EE = 1,14
  TH2F* h_dXvsdY_FH_Acorr[2];  //FH = 1,6                                                                                            
  TH2F* h_dXvsdY_AH_Acorr[2]; //AH 1, 20     

  TH1F* h_dX_EE[2];
  TH1F* h_dY_EE[2];
  TH1F* h_dX_FH[2];
  TH1F* h_dY_FH[2];
  TH1F* h_dX_AH[2];
  TH1F* h_dY_AH[2];



  TH2F *h_xVsy[8];
  TH2F *h_xVsy_corr[8];
  TH2F *h_xVsy_EE[8];
  TH2F *h_xVsy_corr_EE[8];
  TH2F *h_xVsy_FH[8];
  TH2F *h_xVsy_corr_FH[8];
  TH2F *h_xVsy_AH[8];
  TH2F *h_xVsy_corr_AH[8];

  
  TH2F *h_xVsy_SS20_EE_l1[8];
  TH2F *h_xVsy_SS20_EE_l14[8];
  TH2F *h_xVsy_SS20_EE_l1_Acorr[8];
  TH2F *h_xVsy_SS20_EE_l14_Acorr[8];

  TH2F *h_xVsy_EE_l1[8];
  TH2F *h_xVsy_EE_l14[8];

  TH2F *h_xVsy_FH_l1[8];
  TH2F *h_xVsy_FH_l6[8];

  TH2F *h_xVsy_AH_l1[8];
  TH2F *h_xVsy_AH_l20[8];
  
  TH2F *h_xVsy_EE_l1_Acorr[8];
  TH2F *h_xVsy_EE_l14_Acorr[8];

  TH2F *h_xVsy_FH_l1_Acorr[8];
  TH2F *h_xVsy_FH_l6_Acorr[8];

  TH2F *h_xVsy_AH_l1_Acorr[8];
  TH2F *h_xVsy_AH_l20_Acorr[8];


  // TH1F* hist_resp_total_0[85];
  // TH1F* hist_resp_SS_EE_0[85];
  // TH1F* hist_resp_SS_FH_0[85];
  TH1F* hist_resp_total_1[85];
  TH1F* hist_resp_SS_EE_1[85];
  TH1F* hist_resp_SS_FH_1[85];
  TH1F* hist_resp_total_funct[85];
  TH1F* hist_resp_SS_EE_funct[85];
  TH1F* hist_resp_SS_FH_funct[85];
  TH1F* hist_resp_total_trimAhcal[85];
  TH1F* hist_resp_SS_EE_trimAhcal[85];
  TH1F* hist_resp_SS_FH_trimAhcal[85];


  TH1F* hist_total;
  TH1F* hist_frac_EE;
  TH1F* hist_frac_FH;
  TH1F* hist_frac_AH;
  
  
};
#endif

#ifdef AnalyzeHGCOctTB_cxx

void AnalyzeHGCOctTB::BookHistogram(const char *outFileName, const char* conf,  const char* chi_method) {

  conf_ = conf;
  oFile = new TFile(outFileName, "recreate");

  int en = atoi(chi_method);
  float xlow = 0.0;
  float xhigh = en*100*1.5;
  //int xbin = (xhigh-xlow)/6;
  int xbin = 200;//*0.5*en;
  ////////// Book histograms here  //////////////
  int En_bin[8]={20,50,80,100,120,200,250,300};
  int Xmax[3]={10,10,40};
  int Xmin[3]={-10,-10,-40};
  int Xbin[3]={20,20,80};
  char hname[1000];
  const char* compart[3]= {"EE","FH","AH"};
  //histograms for muons
  hist_total = new TH1F("hist_total","total rechits",1000,0,2000);
  hist_frac_EE = new TH1F("hist_frac_EE","",500, 0,500);
  hist_frac_FH=new TH1F("hist_frac_FH","",50, 0,50);
  hist_frac_AH=new TH1F("hist_frac_AH","",50, 0,50);


  for(int ibin=0;ibin<3;ibin++)
    {
      
      sprintf(hname,"h_xVs_trackX_%s",compart[ibin]);
      h_xVs_trackX[ibin]= new TH2F(hname, hname, 20,-10,10,Xbin[ibin],Xmin[ibin],Xmax[ibin]);//,10,-10,10);//Xmin[ibin],Xmax[ibin],80,Xmin[ibin],Xmax[ibin]);
      sprintf(hname,"h_xVs_1trackX_%s",compart[ibin]);
      h_xVs_1trackX[ibin]= new TH2F(hname, hname,20,-10,10, Xbin[ibin],Xmin[ibin],Xmax[ibin]);//,10,-10,10);//80,Xmin[ibin],Xmax[ibin],80,Xmin[ibin],Xmax[ibin]);
      sprintf(hname,"h_yVs_trackX_%s",compart[ibin]);
      h_yVs_trackY[ibin]= new TH2F(hname, hname, 20,-10,10,Xbin[ibin],Xmin[ibin],Xmax[ibin]);//,10,-10,10);//80,Xmin[ibin],Xmax[ibin],80,Xmin[ibin],Xmax[ibin]);
      sprintf(hname,"h_yVs_1trackX_%s",compart[ibin]);
      h_yVs_1trackY[ibin]= new TH2F(hname, hname,20,-10,10, Xbin[ibin],Xmin[ibin],Xmax[ibin]);//,10,-10,10);//80,Xmin[ibin],Xmax[ibin],80,Xmin[ibin],Xmax[ibin]);          
    }
  const char* lay_EE[2]={"layer1","layer14"};
  const	char* lay_FH[2]={"layer1","layer6"};
  const	char* lay_AH[2]={"layer1","layer20"};
  int xmax[3]={20,20,20};
  int xmin[3]={-20,-20,-20};

  for(int ibin=0;ibin<2;ibin++)
    {
      
      sprintf(hname,"h_dXvsdY_%s_%s",compart[0],lay_EE[ibin]);
      h_dXvsdY_EE[ibin]= new TH2F(hname, hname, 80,xmin[ibin],xmax[ibin],80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dXvsdY_%s_%s",compart[1],lay_FH[ibin]);
      h_dXvsdY_FH[ibin]= new TH2F(hname, hname, 80,xmin[ibin],xmax[ibin],80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dXvsdY_%s_%s",compart[2],lay_AH[ibin]);
      h_dXvsdY_AH[ibin]= new TH2F(hname, hname, 80,xmin[ibin],xmax[ibin],80,xmin[ibin],xmax[ibin]);

      sprintf(hname,"h_dXvsdY_Acorr_%s_%s",compart[0],lay_EE[ibin]);
      h_dXvsdY_EE_Acorr[ibin]= new TH2F(hname, hname, 80,xmin[ibin],xmax[ibin],80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dXvsdY_Acorr_%s_%s",compart[1],lay_FH[ibin]);
      h_dXvsdY_FH_Acorr[ibin]= new TH2F(hname, hname, 80,xmin[ibin],xmax[ibin],80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dXvsdY_Acorr_%s_%s",compart[2],lay_AH[ibin]);
      h_dXvsdY_AH_Acorr[ibin]= new TH2F(hname, hname, 80,xmin[ibin],xmax[ibin],80,xmin[ibin],xmax[ibin]);

      sprintf(hname,"h_dX_%s_%s",compart[0],lay_EE[ibin]);
      h_dX_EE[ibin] = new TH1F(hname, hname, 80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dY_%s_%s",compart[0],lay_EE[ibin]);
      h_dY_EE[ibin] = new TH1F(hname, hname, 80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dX_%s_%s",compart[1],lay_FH[ibin]);
      h_dX_FH[ibin] = new TH1F(hname, hname, 80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dY_%s_%s",compart[1],lay_FH[ibin]);
      h_dY_FH[ibin] = new TH1F(hname, hname, 80,xmin[ibin],xmax[ibin]);

      sprintf(hname,"h_dX_%s_%s",compart[2],lay_AH[ibin]);
      h_dX_AH[ibin] = new TH1F(hname, hname, 80,xmin[ibin],xmax[ibin]);
      sprintf(hname,"h_dY_%s_%s",compart[2],lay_AH[ibin]);
      h_dY_AH[ibin] = new TH1F(hname, hname, 80,xmin[ibin],xmax[ibin]);


			       
    }
  
  
  for(int ibin=0;ibin<8;ibin++)
    {
      sprintf(hname,"h_xVsy_beforeCorr_%d",En_bin[ibin]);
      h_xVsy[ibin] = new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_AfterCorr_%d",En_bin[ibin]);
      h_xVsy_corr[ibin]=new TH2F(hname, hname, 80,-40,40,80,-40,40);

      sprintf(hname,"h_xVsy_beforeCorr_EE_%d",En_bin[ibin]);
      h_xVsy_EE[ibin] = new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_AfterCorr_EE_%d",En_bin[ibin]);
      h_xVsy_corr_EE[ibin]=new TH2F(hname, hname, 80,-40,40,80,-40,40);

       sprintf(hname,"h_xVsy_beforeCorr_FH_%d",En_bin[ibin]);
      h_xVsy_FH[ibin] = new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_AfterCorr_FH_%d",En_bin[ibin]);
      h_xVsy_corr_FH[ibin]=new TH2F(hname, hname, 80,-40,40,80,-40,40);
       sprintf(hname,"h_xVsy_beforeCorr_AH_%d",En_bin[ibin]);
      h_xVsy_AH[ibin] = new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_AfterCorr_AH_%d",En_bin[ibin]);
      h_xVsy_corr_AH[ibin]=new TH2F(hname, hname, 80,-40,40,80,-40,40);
      
      sprintf(hname,"h_xVsy_SS20_EE_l1_%d",En_bin[ibin]);
      h_xVsy_SS20_EE_l1[ibin]= new TH2F(hname, hname, 80,-10,10,80,-10,10);
      sprintf(hname,"h_xVsy_SS20_EE_l14_%d",En_bin[ibin]);
      h_xVsy_SS20_EE_l14[ibin]= new TH2F(hname, hname, 80,-10,10,80,-10,10);
      sprintf(hname,"h_xVsy_SS20_EE_l1_Acorr_%d",En_bin[ibin]);
      h_xVsy_SS20_EE_l1_Acorr[ibin]= new TH2F(hname, hname, 80,-10,10,80,-10,10);
      sprintf(hname,"h_xVsy_SS20_EE_l14_Acorr_%d",En_bin[ibin]);
      h_xVsy_SS20_EE_l14_Acorr[ibin]= new TH2F(hname, hname, 80,-10,10,80,-10,10);


      sprintf(hname,"h_xVsy_EE_l1_%d",En_bin[ibin]);
      h_xVsy_EE_l1[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_EE_l14_%d",En_bin[ibin]);
      h_xVsy_EE_l14[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_EE_l1_Acorr_%d",En_bin[ibin]);
      h_xVsy_EE_l1_Acorr[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_EE_l14_Acorr_%d",En_bin[ibin]);
      h_xVsy_EE_l14_Acorr[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);

       sprintf(hname,"h_xVsy_FH_l1_%d",En_bin[ibin]);
      h_xVsy_FH_l1[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_FH_l6_%d",En_bin[ibin]);
      h_xVsy_FH_l6[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_FH_l1_Acorr_%d",En_bin[ibin]);
      h_xVsy_FH_l1_Acorr[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_FH_l6_Acorr_%d",En_bin[ibin]);
      h_xVsy_FH_l6_Acorr[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);

      sprintf(hname,"h_xVsy_AH_l1_%d",En_bin[ibin]);
      h_xVsy_AH_l1[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_AH_l20_%d",En_bin[ibin]);
      h_xVsy_AH_l20[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_AH_l1_Acorr_%d",En_bin[ibin]);
      h_xVsy_AH_l1_Acorr[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);
      sprintf(hname,"h_xVsy_AH_l20_Acorr_%d",En_bin[ibin]);
      h_xVsy_AH_l20_Acorr[ibin]= new TH2F(hname, hname, 80,-40,40,80,-40,40);

      


    }
  h_muo_dxvsdy= new TH2F("h_muo_dxvsdy","dX vs dY",200,-20,20,200,-20,20);
  h_dR = new TH1F("h_dR","sqrt(dx,dy)",500,0,10);

 h_muo_dxvsdy_afterCorr= new TH2F("h_muo_dxvsdy_afterCorr","dX vs dY after correc",200,-20,20,200,-20,20);
 h_dR_afterCorr = new TH1F("h_dR_afterCorr","sqrt(dx,dy) after correc",500,0,10);

   const char *baseline[6]= {"Nocut","FH-noisy","Good-track","MuonVeto","pres-shower","track-impact"};//,"dPhi_Met","Met_100","Met_250","Met_600","st_300_Met100","pt_st_Met_250","st_300_Met250","nocut"};
  char name[100],title[100];
  char hname2[1000],hname1[1000], hname1_2d[1000];
  for(int i=0;i<6;i++)
    {
      cout<<baseline[i]<<endl;
      sprintf(hname,"h_truebeamenergy_after_%s",baseline[i]);
      h_true_beamenergy[i]=new TH1F(hname,hname,400,0,400);
    }
int Elist[85] =  {10, 14 , 18 , 22 , 26 , 30,  34 , 38 , 42,  46,  50,  54,  58,  62,  66,  70,  74,  78,
                     82,  86,  90,  94,  98, 102, 106, 110, 114 ,118, 122, 126, 130, 134, 138, 142 ,146, 150,
                     154 ,158 ,162 ,166 ,170 ,174 ,178, 182, 186, 190 ,194 ,198, 202, 206, 210, 214 ,218, 222,
                     226, 230, 234 ,238, 242 ,246 ,250, 254, 258 ,262 ,266, 270, 274, 278, 282, 286, 290, 294,
                     298, 302 ,306, 310, 314, 318 ,322, 326, 330, 334 ,338 ,342 ,346};

 chi2_mehtod1 =oFile->mkdir("chi2_mehtod0");
 //chi2_mehtod1 =oFile->mkdir("chi2_mehtod1");
   for(int i = 0; i < 85; i++) {
     char* temp = new char[100];
     double xmin = -1;
     double xmax = 10*Elist[i];
     double xbin = 700;
     
     // chi2_mehtod0->cd();
     // sprintf(temp,"Sim_chi2_method0_TrueEn_%d_to_%d",Elist[i],Elist[i+1]);
     // hist_resp_total_0[i]= new TH1F(temp,temp,xbin,0,xmax);
     // sprintf(temp,"Sim_chi2_method0_TrueEn_%d_to_%d_SS_EE",Elist[i],Elist[i+1]);
     // hist_resp_SS_EE_0[i] = new TH1F(temp,temp,xbin,0,xmax);
     // sprintf(temp,"Sim_chi2_method0_TrueEn_%d_to_%d_SS_FH",Elist[i],Elist[i+1]);
     // hist_resp_SS_FH_0[i] = new TH1F(temp,temp,xbin,0,xmax);

     chi2_mehtod1->cd();

     sprintf(temp,"Sim_chi2_method1_TrueEn_%d_to_%d",Elist[i],Elist[i+1]);
     hist_resp_total_1[i]= new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_chi2_method1_TrueEn_%d_to_%d_SS_EE",Elist[i],Elist[i+1]);
     hist_resp_SS_EE_1[i] = new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_chi2_method1_TrueEn_%d_to_%d_SS_FH",Elist[i],Elist[i+1]);
     hist_resp_SS_FH_1[i] = new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_fixwt_method1_TrueEn_%d_to_%d",Elist[i],Elist[i+1]);
     hist_fixwt_total[i]= new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_fixwt_method1_TrueEn_%d_to_%d_SS_EE",Elist[i],Elist[i+1]);
     hist_fixwt_SS_EE[i] = new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_fixwt_method1_TrueEn_%d_to_%d_SS_FH",Elist[i],Elist[i+1]);
     hist_fixwt_SS_FH[i] = new TH1F(temp,temp,xbin,0,xmax);

      sprintf(temp,"Sim_chi2_method1_TrueEn_funct_%d_to_%d",Elist[i],Elist[i+1]);
     hist_resp_total_funct[i]= new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_chi2_method1_TrueEn_funct_%d_to_%d_SS_EE",Elist[i],Elist[i+1]);
     hist_resp_SS_EE_funct[i] = new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_chi2_method1_TrueEn_funct_%d_to_%d_SS_FH",Elist[i],Elist[i+1]);
     hist_resp_SS_FH_funct[i] = new TH1F(temp,temp,xbin,0,xmax);

        sprintf(temp,"Sim_chi2_method1_TrueEn_trimAhcal_%d_to_%d",Elist[i],Elist[i+1]);
     hist_resp_total_trimAhcal[i]= new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_chi2_method1_TrueEn_trimAhcal_%d_to_%d_SS_EE",Elist[i],Elist[i+1]);
     hist_resp_SS_EE_trimAhcal[i] = new TH1F(temp,temp,xbin,0,xmax);
     sprintf(temp,"Sim_chi2_method1_TrueEn_trimAhcal_%d_to_%d_SS_FH",Elist[i],Elist[i+1]);
     hist_resp_SS_FH_trimAhcal[i] = new TH1F(temp,temp,xbin,0,xmax);


     }
   
   //   char* hname = new char[2000];
   for(int i = 0; i < 50; i++) {
     sprintf(hname,"h_weight_%d",i+1);
     
    h_rechit_energy_FB_rel_weightScan[i] = new TH1F(hname,hname,1000,0.0,12000);
    h_rechit_energy_FB_rel_weightScan[i]->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_FB_rel_weightScan[i]->Sumw2();
  }

  h_beamenergy = new TH1F("h_beamenergy","h_beamenergy",320,0,320);
  //  h_true_beamenergy = new TH1F("h_true_beamenergy","h_true_beamenergy",350,0,350);
  h_particle = new TH1F("h_particle","h_particle",600,-300,300);


  h_rechits_mips_EE = new TH1F("h_rechits_mips_EE","h_rechits_mips_EE",500,0,4500);
   h_rechits_mips_FH = new TH1F("h_rechits_mips_FH","h_rechits_mips_FH",500,0,4500);
 h_rechits_mips_AH = new TH1F("h_rechits_mips_AH","h_rechits_mips_AH",500,0,4500);
   
 h_rechits_GeV_EE = new TH1F("h_rechits_GeV_EE","h_rechits_GeV_EE",1000,0,500);
 h_rechits_GeV_FH = new TH1F("h_rechits_GeV_FH","h_rechits_GeV_FH",1000,0,500);
 h_rechits_GeV_AH = new TH1F("h_rechits_GeV_AH","h_rechits_GeV_AH",1000,0,500);
 h_rechitvslambda = new TH2F("h_rechitvslambda","h_rechitvslambda",60,0,12,500,0,3000);

  h_EnergySum_inEE=new TH1F("h_EnergySum_SSinEE","h_EnergySum_SSinEE",500,0,40000);
  h_ssinEE = new TH1F("h_ssinEE","h_ssinEE",300,0,300);
  h_mipsinEE = new TH1F("h_mipsinEE","h_mipsinEE",300,0,300);
  h_EnergySum_inFH=new TH1F("h_EnergySum_SSinFH","h_EnergySum_SSinFH",500,0,40000);
  h_EnergySum_inAH=new TH1F("h_EnergySum_SSinAH","h_EnergySum_SSinAH",500,0,40000);
  h_EnergySum_inEE_inGeV=new TH1F("h_EnergySum_SSinEE_inGeV","h_EnergySum_SSinEE_inGeV",1000,0,1000);
  h_EnergySum_inFH_inGeV=new TH1F("h_EnergySum_SSinFH_inGeV","h_EnergySum_SSinFH_inGeV",1000,0,1000);
  h_EnergySum_inAH_inGeV=new TH1F("h_EnergySum_SSinAH_inGeV","h_EnergySum_SSinAH_inGeV",1000,0,1000);
  h_EnergySum_inGeV = new TH1F("h_EnergySum_inGeV","h_EnergySum_inGeV",1000,0,1000);
  h_EnergySum_ratio_inGeV = new TH1F("h_EnergySum_ratio_inGeV","h_EnergySum_ratio_inGeV",1000,0,10);
  h_EnergySum_ratio_vs_true = new TH2F("h_EnergySum_ratio_vs_true"," ratio vs true",1000,0,10,500,0,500);
  h_EnergySum_ratio_flipped = new TH1F("h_EnergySum_ratio_flipped","pred/true",10000,-10,3);
  h_EnergySum_ratio_vs_true_flipped = new TH2F("h_EnergySum_ratio_vs_true_flipped"," pred/true vs true",1000,0,10,500,0,500);
  h_EnergySum_ratio_inMips = new TH1F("h_EnergySum_ratio_inMips","h_EnergySum_ratio_inMips",10000,0,10);
  h_EnergySum_inGeVvstrue = new TH2F("h_EnergySum_inGeVvstrue","h_EnergySum_inGeV",1000,0,1000,500,0,500);
  h_EnergySum_inMipsvstrue = new TH2F("h_EnergySum_inMipsvstrue","h_EnergySum_inGeV",1000,0,60000,500,0,500);
  h_EnergySum_inAH_vs_true = new TH2F("h_EnergySum_inAH_vs_true","h_EnergySum_inAH_vs_true",800,0,800,500,0,500);
  h_EnergySum_inEE_vs_true = new TH2F("h_EnergySum_inEE_vs_true","h_EnergySum_inEE_vs_true",800,0,800,500,0,500);
  h_EnergySum_inFH_vs_true = new TH2F("h_EnergySum_inFH_vs_true","h_EnergySum_inFH_vs_true",800,0,800,500,0,500);
  h_EnergySum_inAH_Mipsvs_true = new TH2F("h_EnergySum_inAH_Mipsvs_true","h_EnergySum_inAH_Mipsvs_true",1000,0,40000,500,0,500);
  h_EnergySum_inEE_Mipsvs_true = new TH2F("h_EnergySum_inEE_Mipsvs_true","h_EnergySum_inEE_Mipsvs_true",1000,0,60000,500,0,500);
  h_EnergySum_inFH_Mipsvs_true = new TH2F("h_EnergySum_inFH_Mipsvs_true","h_EnergySum_inFH_Mipsvs_true",1000,0,20000,500,0,500);
  h_EnergySum_inEE_vs_FH = new TH2F("h_EnergySum_inEE_vs_FH","enegry in EE vs energy in FH in mips",1000,0,60000,1000,0,20000);
  h_EnergySum_inEE_vs_FH_inGeV = new TH2F("h_EnergySum_inEE_vs_FH_inGeV","enegry in EE vs energy in FH in GeV",1000,0,1000,1000,0,1000);
  h_EnergySum_Rechit = new TH1F("h_EnergySum_Rechit","h_EnergySum_Rechit",500,0,60000);
  h_EnergySum_Rechit_normalized = new TH1F("h_EnergySum_Rechit_normalized","h_EnergySum_Rechit_normalized",500,0,1000);
  h_nrechits = new TH1F("h_nrechits","total rechits in the prototype",1000,0,10000);
  h_nrechits_EE = new TH1F("h_nrechits_EE","rechits in EE",1000,0,10000);
  h_nrechits_FH = new TH1F("h_nrechits_FH","rechits in FH",1000,0,10000);
  h_nrechits_AH = new TH1F("h_nrechits_AH","rechits in AH",1000,0,10000);
  h_ratio_nrechits_AH = new TH1F("h_ratio_nrechits_AH","h_ratio_nrechits_AH",1000,0,10);
  h_ratio_nrechits_FH =	new TH1F("h_ratio_nrechits_FH","h_ratio_nrechits_FH",1000,0,10);
  h_ratio_nrechits_EE =	new TH1F("h_ratio_nrechits_EE","h_ratio_nrechits_EE",1000,0,10);
  h_nrechits_vs_EE = new TH2F("h_nrechits_vs_EE","total vs rechits in EE",1000,0,10000,1000,0,10000);
  h_nrechits_vs_FH = new TH2F("h_nrechits_vs_FH","total vs rechits in FH",1000,0,10000,1000,0,10000);
  h_nrechits_vs_AH = new TH2F("h_nrechits_vs_AH","total vs rechits in AH",1000,0,10000,1000,0,10000);
  h_nrechits_EE_vs_FH = new TH2F("h_nrechits_EE_vs_FH","EE vs rechits in FH",1000,0,10000,1000,0,10000);

}

void AnalyzeHGCOctTB::Alignment_Map_Init() {
  /// Alignment map for config 1, it needs to be changed accorfing to configurations
  
  char* f_name = new char[200];
  sprintf(f_name,"../../Alignment_Map.txt");
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
  sprintf(f_name,"../../Noise_Map.txt");
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
    sprintf(f_name,"../../config_maps/moduleMAP_config1.txt");
    cout<<"\n\nINFO: Mapping module configuration ALPHA (oct10-oct17) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else if(strcmp(config,"bravo")==0 || strcmp(config,"config2")==0) {
    sprintf(f_name,"../../config_maps/moduleMAP_config2.txt");
    cout<<"\n\nINFO: Mapping module configuration BRAVO (17oct-22oct) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else if(strcmp(config,"charlie")==0  || strcmp(config,"config3")==0) {
    sprintf(f_name,"../../config_maps/moduleMAP_config3.txt");
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
  sprintf(f_name,"../../config1_lengths.txt");

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



void AnalyzeHGCOctTB::Chi2_Weight_Map_Init(int chi2_method) {
  char *f_name_EH = new char[2000];
  char *f_name_H = new char[2000];
  char *method = new char[2000];
  bool UseSimWeights = false;

  switch(chi2_method) {
  case 0 : sprintf(method,"%d : rechit energy sum as input to chi2",chi2_method);  break;
  case 1 : sprintf(method,"%d : detector scale with No offset for H hadrons as input to chi2",chi2_method);  break;
  case 2 : sprintf(method,"%d : detector scale with offset (0.4 GeV) for H hadrons as input to chi2",chi2_method);  break;
  case 3 : sprintf(method,"%d : detector scale with offset (0.4 GeV) for H hadrons AND events around the core of beam energy as input to chi2",chi2_method);  break;
  default : cout<<"Incorrect method input use between 0-3"<<endl; exit(0);
}
  sprintf(f_name_EH,"./txt_maps/chi2_flatEn/sim/chi2_calibFact_EH_hadrons_flatEn_%d_scalMC_2sigma.txt",chi2_method);
  sprintf(f_name_H,"./txt_maps/chi2_flatEn/sim/chi2_calibFact_H_hadrons_flatEn_%d_scalMC_2sigma.txt",chi2_method);



  std::ifstream in_EH(f_name_EH);
  std::ifstream in_H(f_name_H);
  if(!in_EH) {
    cout<<"ERROR => "<<f_name_EH<<" Not found"<<endl;
    //return;                                                                                                                                
    exit(0);
  }
  if(!in_H) {
    cout<<"ERROR => "<<f_name_H<<" Not found"<<endl;
    //return;                                                                                                                                
    exit(0);
  }
int beamEnergy;
 float  w1, w2, w3,ew1,ew2,ew3;

  cout<<"File name = "<<f_name_EH<<endl;
  while(in_EH>>beamEnergy>>w1>>w2>>w3>>ew1>>ew2>>ew3){
    std::pair<int, std::vector<float>> temp_pair;
    std::vector<float> temp_vector;
    temp_vector.push_back(w1);
    temp_vector.push_back(w2);
    temp_vector.push_back(w3);


    temp_pair = std::make_pair(beamEnergy,temp_vector);
    chi2_weights_EH.insert(temp_pair);
  }
  cout<<BOLDGREEN<<"INFO: Chi2 calibration Map initialized for EH hadrons!! "<<RESET<<endl;
  beamEnergy = -1.0; w1 = -1.0; w2 = -1.0; w3 = -1.0;


  cout<<"File name = "<<f_name_H<<endl;
  while(in_H>>beamEnergy>>w1>>w2>>w3>>ew1>>ew2>>ew3){
    std::pair<int, std::vector<float>> temp_pair;
    std::vector<float> temp_vector;
    temp_vector.push_back(w1);
    temp_vector.push_back(w2);
    temp_vector.push_back(w3);


    temp_pair = std::make_pair(beamEnergy,temp_vector);
    chi2_weights_H.insert(temp_pair);
  }
   cout<<BOLDGREEN<<"INFO: Chi2 calibration Map initialized for H hadrons!! "<<RESET<<endl<<endl;

  cout<<BOLDGREEN<<"INFO: Chi2 minimization method used =>  "<<method<<RESET<<endl<<endl;

  cout<<chi2_weights_H.size()<<"\t"<<chi2_weights_EH.size()<<endl;
}
AnalyzeHGCOctTB::AnalyzeHGCOctTB(const TString &inputFileList, const char *outFileName, const char* dataset, const char* config, const char*chi_method,const char* folderID) {
  TChain *tree = new TChain("pion_variables");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  int chi2_method= atoi(chi_method);
  HGCNtupleVariables::Init(tree);//, tree2, tree3);                                                                                         
  BookHistogram(outFileName, config, chi_method);
  moduleMap_init(config);
  Alignment_Map_Init();
  Noise_Map_Init();
  layerPosition_init();
  Chi2_Weight_Map_Init(chi2_method);
}


Bool_t AnalyzeHGCOctTB::FillChain(TChain *chain,  const TString &inputFileList) {//TChain *chain2, TChain *chain3, const TString &inputFileList) {

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
    // chain2->Add(buffer.c_str());
    // chain3->Add(buffer.c_str());

  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain3 : " << chain3->GetEntries() << std::endl;

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

  // if (!fChain2) return -5;
  // Long64_t centry2 = fChain2->LoadTree(entry);
  // if (centry2 < 0) return centry2;
  // if (!fChain2->InheritsFrom(TChain::Class()))  return centry2;
  // TChain *chain2 = (TChain*)fChain2;
  // if (chain2->GetTreeNumber() != fCurrent) {
  //   fCurrent = chain->GetTreeNumber();
  //   //    Notify();
  // }

  // if (!fChain3) return -5;
  // Long64_t centry3 = fChain3->LoadTree(entry);
  // if (centry3 < 0) return centry3;
  // if (!fChain3->InheritsFrom(TChain::Class()))  return centry3;
  // TChain *chain3 = (TChain*)fChain3;
  // if (chain3->GetTreeNumber() != fCurrent) {
  //   fCurrent = chain->GetTreeNumber();
  //   //    Notify();
  // }
  
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
