//---------------------------------Spheral++----------------------------------//
// SuperGaussianKernel -- The super gaussian interpolation kernel.
//
// Created by JMO, Wed Dec  1 22:38:21 PST 1999
//----------------------------------------------------------------------------//
#include <math.h>

#include "SuperGaussianKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    // template class SuperGaussianKernel<Dim<1> >;
    // template class SuperGaussianKernel<Dim<2> >;
    // template class SuperGaussianKernel<Dim<3> >;

    template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
    template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
    template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;

    template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;
    template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;
    template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;

//     template<typename Dimension>
//     const double SuperGaussianKernel<Dimension>::mKW = 0.5*(double(Dimension::nDim) + 2.0);

//     template<typename Dimension>
//     const double SuperGaussianKernel<Dimension>::mKGW = 0.5*double(Dimension::nDim);

  }
}

// Define the static volume normalization.
// ASPH SuperGaussian interpolation kernel:
//     1D normalization:   A = Det(H)/pi^(1/2)
//     2D normalization:   A = Det(H)/pi
//     3D normalization:   A = Det(H)/pi^(3/2)
// const double SuperGaussianKernel<Dim<1> >::mVolumeNormalization = 1.0/sqrt(M_PI);
// const double SuperGaussianKernel<Dim<2> >::mVolumeNormalization = 1.0/(M_PI);
// const double SuperGaussianKernel<Dim<3> >::mVolumeNormalization = 1.0/pow(M_PI, 1.5);

// template<typename Dimension>
// const double SuperGaussianKernel<Dimension>::mKernelExtent = 3.0;
