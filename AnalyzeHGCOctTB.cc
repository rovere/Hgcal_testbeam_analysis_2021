#define AnalyzeHGCOctTB_cxx
#include "AnalyzeHGCOctTB.h"

#include <TF1.h>
#include <math.h>

#include <cstring>
#include <iostream>
#include <vector>

#include "CLUEAlgo.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
using namespace std;

// chip 3022,44,3028

int main(int argc, char *argv[])  //, int argvv[])
{
  if (argc < 5) {
    cerr << "Please give 5 arguments "
         << "runList "
         << " "
         << "outputFileName"
         << " "
         << "dataset"
         << " "
         << "configuration"
         << " "
         << "energy" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName = argv[2];
  const char *data = argv[3];
  const char *config = argv[4];
  const char *energy = argv[5];
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy);
  cout << "dataset " << data << " " << endl;
  cout << "configuration " << config << " " << endl;
  cout << "energy  " << energy << " " << endl;
  // int chi2_method = atoi(energy);
  hgcOctTB.EventLoop(data, energy);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data, const char *energy) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  // Long64_t nentries3 = fChain3->GetEntriesFast();
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
  if (!strcmp(conf_, "alpha") || !strcmp(conf_, "config1")) {
    TOTAL_ACTIVE_LAYER = 79;
    EE_LAYER = 28;
    FH_LAYER = 12;
    AH_LAYER = 39;
  } else if (!strcmp(conf_, "bravo") || !strcmp(conf_, "config2")) {
    TOTAL_ACTIVE_LAYER = 78;
    EE_LAYER = 28;
    FH_LAYER = 11;
    AH_LAYER = 39;
  } else if (!strcmp(conf_, "charlie") || !strcmp(conf_, "config3")) {
    TOTAL_ACTIVE_LAYER = 59;
    EE_LAYER = 8;
    FH_LAYER = 12;
    AH_LAYER = 39;

  } else {
    cout << "ERROR: Unknown configuration!!!!" << endl;
    return;
  }

  // counter
  int nHgc = 0, nAhc = 0, nrechits = 0;
  float offset = 159.4;  // A constant to be added to AHC layer z positions
  // cout<<energy<<endl;
  float FH_AH_relative_scale = 0.4;
  float alpha_ = FH_AH_relative_scale;
  float EE_scale = 94.624;     // MIPs per Ge
  float FH_AH_scale = 12.788;  // MIPs per GeV
  float ee_rescaling = 1.035;
  float fh_rescaling = 1.095;
  float ah_rescaling = 1.095;
  if (!strcmp(data, "data")) {
    ee_rescaling = 1;
    fh_rescaling = 1;
    ah_rescaling = 1;
  }
  cout << "sim rescaling"
       << "\t" << ee_rescaling << "\t" << fh_rescaling << "\t" << ah_rescaling
       << endl;

  int Elist[85] = {
      10,  14,  18,  22,  26,  30,  34,  38,  42,  46,  50,  54,  58,  62,  66,
      70,  74,  78,  82,  86,  90,  94,  98,  102, 106, 110, 114, 118, 122, 126,
      130, 134, 138, 142, 146, 150, 154, 158, 162, 166, 170, 174, 178, 182, 186,
      190, 194, 198, 202, 206, 210, 214, 218, 222, 226, 230, 234, 238, 242, 246,
      250, 254, 258, 262, 266, 270, 274, 278, 282, 286, 290, 294, 298, 302, 306,
      310, 314, 318, 322, 326, 330, 334, 338, 342, 346};

  float comb_z_boundaries[50] = {14., 15., 17., 18., 20., 21., 23., 24., 26.,
  27., 29., 30., 32., 33., 35., 36., 37.5, 38.5, 40.5, 41., 43.5, 44.5, 47., 48.,
  50., 51., 53., 54.5, 65., 72., 79., 86., 92., 99., 116., 124., 132., 140., 146., 154.,
  170., 180., 190., 200., 210., 220., 232., 242., 252., 262.};

  // to select the evnts within 2sigma region

  if (DEBUG) cout << "DEBUG: Configuration = " << conf_ << endl;
  if (DEBUG)
    cout << "DEBUG: TOTAL_ACTIVE_LAYER = " << TOTAL_ACTIVE_LAYER << endl;
  if (DEBUG) cout << "DEBUG: EE_LAYER = " << EE_LAYER << endl;
  if (DEBUG) cout << "DEBUG: FH_LAYER = " << FH_LAYER << endl;
  if (DEBUG) cout << "DEBUG: AH_LAYER = " << AH_LAYER << endl;

  if (DEBUG) cout << "DEBUG: Entering event Loop" << endl;

  Long64_t jentry;

  float lambda[79];
  for (int i = 0; i < 79; i++) {
    lambda[i] = layer_positions[i + 1].at(
        2);  // for nuclear interaction length  & //pi_lambda use -at(3)
    //     cout<<lambda[i]<<endl;
  }
  int ahcal_layer[10] = {43, 47, 51, 55, 59, 63,
                         67, 71, 75, 79};  // selected layers (10 in total) out
                                           // of AHCAL (39 in total)

  for (jentry = 0; jentry < nentries; jentry++, hgc_jentry++) {
    //   for (jentry=0; jentry<10000;jentry++,hgc_jentry++) {
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int(progress);
    if (k > decade) cout << 10 * k << " %" << endl;
    decade = k;

    // ===============read this entry == == == == == == == == == == ==

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) {
      break;
      cout << "Breaking" << endl;
    }
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // cout<<"insides the event loop: check1"<<endl;

    ////   MESSAGE ////
    // apparently code is not running beyond this point ///
    event_count[0]++;

    event_count[1]++;
    h_true_beamenergy[4]->Fill(trueBeamEnergy);
    totalEnergy_inGeV = 0;
    EnergySum_SSinEE = 0;
    EnergySum_SSinFH = 0;
    Esum_rechits_FH = 0;
    Esum_rechits_EE = 0;
    Esum_rechits_AH = 0;
    Esum_rechits_FH_inGeV = 0;
    Esum_rechits_EE_inGeV = 0;
    Esum_rechits_AH_inGeV = 0;

    ////////////////////////////////////////////
    //            HGCAL Part                  //
    ////////////////////////////////////////////
    total_rechits = 0;
    rechits_EE = 0;
    rechits_FH = 0;
    rechits_AH = 0;
    // if(event ==3482)
    // {//continue;
    int nrechit_trimAhcal = 0;
    /// Read HGCAL + AHCAL (only 10 layers out of 39) combined Tree
    if (DEBUG) cout << "DEBUG: Start Analylizing HGCAL  RecHits!!" << endl;
    float rechitEnergySum_AH = 0.0;
    for (int i = 0; i < rechitEn_trimAhcal->size(); i++) {
      float energy = rechitEn_trimAhcal->at(i);
      if (comb_rechit_z_trimAhcal->at(i) < 54)  // selecting EE rechits
      {
        Esum_rechits_EE += (energy / ee_rescaling);
        rechits_EE++;
      } else if ((comb_rechit_z_trimAhcal->at(i) > 54) &&
                 (comb_rechit_z_trimAhcal->at(i) <
                  154))  // Selecting FH rechits
      {
        Esum_rechits_FH += (energy / fh_rescaling);
      } else if (comb_rechit_z_trimAhcal->at(i) > 154)  // selecting AH rechits
      {
        rechitEnergySum_AH += (energy / ah_rescaling);
      }
    }  // Nrechits loop

    h_true_beamenergy[5]->Fill(trueBeamEnergy);
    nrechit_trimAhcal = 0;
    rechits_AH = nAhc;
    total_rechits = nrechits;
    // convert the energy into GeV using detector level calibration
    Esum_rechits_EE_inGeV = 0.0105 * Esum_rechits_EE;
    Esum_rechits_FH_inGeV = 0.0812 * Esum_rechits_FH;  // updated relative
                                                       // weight
    Esum_rechits_AH_inGeV = 0.12508 * rechitEnergySum_AH;
    auto Esum_allRecHits_inGeV =
        Esum_rechits_AH_inGeV + Esum_rechits_FH_inGeV + Esum_rechits_EE_inGeV;
    if (DEBUG) cout << "DEBUG: End of Entry = " << jentry << endl;
    if (DEBUG) cout << "DEBUG: ****************** " << endl;
    //      cout<<"\t"<<"beforE"<<endl;

    // Compute clusters using CLUE
    std::array<LayerTiles, NLAYERS> tiles;
    PointsCloud pcloud;
    constexpr float dc = 1.3;
    constexpr float rhoc = 2.f;
    constexpr float outlierDeltaFactor = 2.f;
    pcloud.x = *comb_rechit_x_trimAhcal;
    pcloud.y = *comb_rechit_y_trimAhcal;
    pcloud.z = *comb_rechit_z_trimAhcal;
    pcloud.weight = *rechitEn_trimAhcal;
    updateLayers(pcloud, comb_z_boundaries);
//    pcloud.layer = *rechit_layer;
    pcloud.outResize(rechit_x->size());

    compute_histogram(tiles, pcloud);
    calculate_density(tiles, pcloud, dc);
    calculate_distanceToHigher(tiles, pcloud, outlierDeltaFactor, dc);
    auto total_clusters = findAndAssign_clusters(pcloud, outlierDeltaFactor, dc, 2.);
    auto clusters = getClusters(total_clusters, pcloud);
    float total_energy_clustered = 0.f;
    for (auto const & cl : clusters) {
      auto pos = cl.position();
      std::cout << std::get<0>(pos) << " "
        << std::get<1>(pos) << " "
        << std::get<2>(pos) << " "
        << std::endl;
      total_energy_clustered += cl.energy();
    }
    if (1) {
      std::cout << "True energy: " << trueBeamEnergy
                << " Sum(RecHits_energy): " << Esum_allRecHits_inGeV
                << " Ratio: " << Esum_allRecHits_inGeV / trueBeamEnergy
                << " Clustered Energy: " << total_energy_clustered
                << " Ratio: " << total_energy_clustered / Esum_allRecHits_inGeV
                << " Total number of clusters: " << total_clusters
                << std::endl;
    }
  }  // loop over entries

  ///////////////////////////////////////////////////////////
  ///////  E N D     O F     E N T R Y     L O O P     //////
  ///////////////////////////////////////////////////////////
  //  gSystem->Exit(0);
  cout << "Got Out " << jentry << endl;
}
