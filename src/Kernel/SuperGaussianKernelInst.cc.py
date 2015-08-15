text = """
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
"""

if ndim == 1:
    text += """
    template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
    template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;
"""

elif ndim == 2:
    text += """
    template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
    template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;
"""

elif ndim == 3:
    text += """
    template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;
    template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;
"""

text += """
    template class SuperGaussianKernel<Dim< %(ndim)s > >;
  }
}
"""
