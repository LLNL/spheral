text = """
//---------------------------------Spheral++----------------------------------//
// QuinticSplineKernel -- A quintic spline, as described in
// Bonet & Kulasegaruam 2002, Appl. Math. Comput., 126, 133-155.
//
// Kernel extent: 2.0
//
// Created by JMO, Wed Jul  9 16:24:25 PDT 2008
//----------------------------------------------------------------------------//
#include "QuinticSplineKernel.cc"


namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
"""

if ndim == "1":
   text += """
template<>
QuinticSplineKernel< Dim<1> >::QuinticSplineKernel():
  Kernel<Dim<1>, QuinticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(1.0/16.0);
  setKernelExtent(2.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}
"""

elif ndim == "2":
    text += """
template<>
QuinticSplineKernel< Dim<2> >::QuinticSplineKernel():
  Kernel<Dim<2>, QuinticSplineKernel< Dim<2> > >() {
  setVolumeNormalization(3.0/(16.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}
"""

else:
    assert ndim == "3"
    text += """
template<>
QuinticSplineKernel< Dim<3> >::QuinticSplineKernel():
  Kernel<Dim<3>, QuinticSplineKernel< Dim<3> > >() {
  setVolumeNormalization(7.0/(40.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}
"""

text += """

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class QuinticSplineKernel< Dim< %(ndim)s >  >;
}
}
"""
