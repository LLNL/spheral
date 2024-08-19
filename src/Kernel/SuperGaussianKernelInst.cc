//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include <math.h>

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
"""

if ndim == "1":
    text += """
    template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
    template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;
"""

elif ndim == "2":
    text += """
    template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
    template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;
"""

elif ndim == "3":
    text += """
    template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;
    template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;
"""

text += """
    template class SuperGaussianKernel<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
"""

if ndim == "1":
    text += """
    template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
    template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;
"""

elif ndim == "2":
    text += """
    template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
    template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;
"""

elif ndim == "3":
    text += """
    template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;
    template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;
"""

text += """
    template class SuperGaussianKernel<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
"""

if ndim == "1":
    text += """
    template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
    template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;
"""

elif ndim == "2":
    text += """
    template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
    template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;
"""

elif ndim == "3":
    text += """
    template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;
    template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;
"""

text += """
    template class SuperGaussianKernel<Dim<3> >;
#endif
}