//---------------------------------Spheral++----------------------------------//
// KernelTraits -- trait class for dimensional type information.
//
// Created by JMO, Mon Nov 29 21:51:58 PST 1999
//----------------------------------------------------------------------------//

#ifndef KERNELTRAITS_HH
#define KERNELTRAITS_HH

#include "Geometry/Dimension.hh"
#include "BSplineKernel.hh"

template<typename KernelType> 
class KernelTraits {
};

// ASPH Spline interpolation kernel:
//    1D normalization:   A = 2*Det(H)/3
//    2D normalization:   A = 10*Det(H)/(7*pi)
//    3D normalization:   A = Det(H)/pi

template<>
class KernelTraits<BSplineKernel<Dimension<1> > > {
public:
  const double kernelNormalization() const {
    return 0.666666666667;
  }
};

template<>
class KernelTraits<BSplineKernel<Dimension<2> > > {
public:
  const double kernelNormalization() const {
    return 0.454728408834;
  }
};

template<>
class KernelTraits<BSplineKernel<Dimension<3> > > {
public:
  const double kernelNormalization() const {
    return 0.318309886184;
  }
};

#endif
