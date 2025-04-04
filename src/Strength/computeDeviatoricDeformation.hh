//---------------------------------Spheral++----------------------------------//
// computeDeviatoricDeformation
//
// Function to generalize computing the deviatoric deformation across solid
// hydrodyanmics packages taking into account how to fill in the missing
// strain components in lower dimensions.
//
// Compute the deviatoric deformation, where for dimensions < 3 we either use
// plane-strain : all non-evolved strains are zero   (planeStrain = true)
// plane-stress : all non-evolved stresses are zero  (planeStrain = false)
//
// Created by JMO, Thu Apr  3 14:41:47 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_computeDeviatoricDeformation__
#define __Spheral_computeDeviatoricDeformation__

namespace Spheral {

template<typename SymTensor>
inline
SymTensor
computeDeviatoricDeformation(const SymTensor& deformation,
                             const bool planeStrain) {
  return deformation - deformation.Trace()/3.0*SymTensor::one;
  // const auto Dtrace = (planeStrain ?
  //                      deformation.Trace() :
  //                      0.0);
  // return deformation - Dtrace/3.0*SymTensor::one;
}

// // In 3D the plane-strain vs. plane-stress distinction is not used
// template<>
// inline
// GeomSymmetricTensor<3>
// computeDeviatoricDeformation<GeomSymmetricTensor<3>>(const GeomSymmetricTensor<3>& deformation,
//                                                      const bool planeStrain) {
//   return deformation - deformation.Trace()/3.0*GeomSymmetricTensor<3>::one;
// }

}

#endif
