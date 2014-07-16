//---------------------------------Spheral++----------------------------------//
// GaussianKernel -- The gaussian interpolation kernel.
//
// Created by JMO, Wed Dec  1 14:38:51 PST 1999
//----------------------------------------------------------------------------//
#include <math.h>

#include "Kernel.hh"
#include "GaussianKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class GaussianKernel<Dim<1> >;
    template class GaussianKernel<Dim<2> >;
    template class GaussianKernel<Dim<3> >;
  }
}


// Define the static volume normalization.
// ASPH Gaussian interpolation kernel:
//     1D normalization:   A = Det(H)/sqrt(pi)
//     2D normalization:   A = Det(H)/pi
//     3D normalization:   A = Det(H)/pi**1.5
// const double GaussianKernel<Dim<1> >::mVolumeNormalization = 1.0/sqrt(M_PI);
// const double GaussianKernel<Dim<2> >::mVolumeNormalization = 1.0/M_PI;
// const double GaussianKernel<Dim<3> >::mVolumeNormalization = 1.0/pow(M_PI, 1.5);

// template<typename Dimension>
// const double GaussianKernel<Dimension>::mKernelExtent = 3.0;
