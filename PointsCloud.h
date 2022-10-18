#ifndef Points_Cloud_h
#define Points_Cloud_h

#include <vector>

struct PointsCloud {
  PointsCloud() = default;

  void outResize(unsigned int const& nPoints) {
    rho.resize(nPoints);
    delta.resize(nPoints);
    nearestHigher.resize(nPoints);
    clusterIndex.resize(nPoints);
    followers.resize(nPoints);
    isSeed.resize(nPoints);
    n = nPoints;
  }

  // Input variables
  std::vector<float> x;
  std::vector<float> y;
  std::vector<unsigned int> layer;
  std::vector<float> weight;

  // Output variables
  std::vector<float> rho;
  std::vector<float> delta;
  std::vector<int> nearestHigher;
  std::vector<std::vector<int>> followers;
  std::vector<int> isSeed;
  std::vector<int> clusterIndex;

  unsigned int n;
};

#endif
