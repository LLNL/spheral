//---------------------------------Spheral++----------------------------------//
// W4SplineKernel -- The B spline interpolation kernel.
//
// Created by JMO, Mon Nov 29 22:57:26 PST 1999
//----------------------------------------------------------------------------//

#include <math.h>

#include "Kernel.hh"
#include "W4SplineKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class W4SplineKernel<Dim<1> >;
    template class W4SplineKernel<Dim<2> >;
    template class W4SplineKernel<Dim<3> >;
  }
}

// Define the static volume normalization.
// ASPH Spline interpolation kernel:
//    1D normalization:   A = 1*Det(H)
//    2D normalization:   A = 30*Det(H)/(7*pi)
//    3D normalization:   A = 5*Det(H)/(6*pi)
// const double W4SplineKernel<Dim<1> >::mVolumeNormalization = 1.0;
// const double W4SplineKernel<Dim<2> >::mVolumeNormalization = 30.0/(7.0*M_PI);
// const double W4SplineKernel<Dim<3> >::mVolumeNormalization = 5.0/(6.0*M_PI);

// template<typename Dimension>
// const double W4SplineKernel<Dimension>::mKernelExtent = 2.0;
