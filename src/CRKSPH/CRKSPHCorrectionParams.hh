//Parameters used for CRK kernels. 
#ifndef __Spheral_CRKSPHCorrectionParams_hh__
#define __Spheral_CRKSPHCorrectionParams_hh__
//Enumerated type for the corrected Kernels
namespace Spheral {
namespace CRKSPHSpace {

enum CRKOrder {//Used to assign the order of the corrections
  ZerothOrder = 0,
  LinearOrder = 1,
  QuadraticOrder = 2
};

enum CRKVolumeType { // Choices for the CRK volume weighting
  CRKMassOverDensity = 0,
  CRKSumVolume = 1,
  CRKVoronoiVolume = 2,
  CRKHullVolume = 3,
  HVolume = 4,
};

}
}
#endif
