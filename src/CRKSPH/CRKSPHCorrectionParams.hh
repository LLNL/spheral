//Parameters used for CRK kernels. 
#ifndef __Spheral_CRKSPHCorrectionParams_hh__
#define __Spheral_CRKSPHCorrectionParams_hh__
//Enumerated type for the corrected Kernels
namespace Spheral {

enum class CRKOrder : int {//Used to assign the order of the corrections
  ZerothOrder = 0,
  LinearOrder = 1,
  QuadraticOrder = 2,
  CubicOrder = 3,
  QuarticOrder = 4,
  QuinticOrder = 5,
  SexticOrder = 6,
  SepticOrder = 7,
};

enum class CRKVolumeType : int { // Choices for the CRK volume weighting
  CRKMassOverDensity = 0,
  CRKSumVolume = 1,
  CRKVoronoiVolume = 2,
  CRKHullVolume = 3,
  HVolume = 4,
};

}
#endif
