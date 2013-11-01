//---------------------------------Spheral++----------------------------------//
// PiGaussianKernel -- The gaussian interpolation kernel.
//
// Created by JMO, Wed Dec  1 14:38:51 PST 1999
//----------------------------------------------------------------------------//
#include <math.h>

#include "Kernel.hh"
#include "PiGaussianKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class PiGaussianKernel<Dim<1> >;
    template class PiGaussianKernel<Dim<2> >;
    template class PiGaussianKernel<Dim<3> >;
  }
}

// Define the static volume normalization.
// ASPH Pi Gaussian interpolation kernel:
//     1D normalization:   A = 2*K^(1/4)*Det(H)/Gamma(1/4)
//     2D normalization:   A = 2*K^(1/2)*Det(H)/pi^(3/2)
//     3D normalization:   A = K^(3/4)*Det(H)/(pi*Gamma(3/4))
// const double PiGaussianKernel<Dim<1> >::mVolumeNormalization = 2.0/3.6256099082;
// const double PiGaussianKernel<Dim<2> >::mVolumeNormalization = 2.0/pow(M_PI, 1.5);
// const double PiGaussianKernel<Dim<3> >::mVolumeNormalization = 1.0/(M_PI*1.2254167024);

// template<typename Dimension>
// const double PiGaussianKernel<Dimension>::mKernelExtent = 2.0;
