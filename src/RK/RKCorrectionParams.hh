//Parameters used for RK kernels. 
#ifndef __Spheral_RKCorrectionParams_hh__
#define __Spheral_RKCorrectionParams_hh__
//Enumerated type for the corrected Kernels
namespace Spheral {

enum class RKOrder : int {//Used to assign the order of the corrections
  ZerothOrder = 0,
  LinearOrder = 1,
  QuadraticOrder = 2,
  CubicOrder = 3,
  QuarticOrder = 4,
  QuinticOrder = 5,
  SexticOrder = 6,
  SepticOrder = 7,
};

enum class RKVolumeType : int { // Choices for the RK volume weighting
  RKMassOverDensity = 0,
  RKSumVolume = 1,
  RKVoronoiVolume = 2,
  RKHullVolume = 3,
  HVolume = 4,
};

}
#endif
