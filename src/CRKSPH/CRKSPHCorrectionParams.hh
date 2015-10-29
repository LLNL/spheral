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
}
}
#endif
