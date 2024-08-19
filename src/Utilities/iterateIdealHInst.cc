//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/iterateIdealH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template void iterateIdealH<Dim<1> >(DataBase<Dim<1> >&,
                                       const vector<Boundary<Dim<1> >*>&,
                                       const TableKernel<Dim<1> >&,
                                       const SmoothingScaleBase<Dim<1> >&,
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template void iterateIdealH<Dim<2> >(DataBase<Dim<2> >&,
                                       const vector<Boundary<Dim<2> >*>&,
                                       const TableKernel<Dim<2> >&,
                                       const SmoothingScaleBase<Dim<2> >&,
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template void iterateIdealH<Dim<3> >(DataBase<Dim<3> >&,
                                       const vector<Boundary<Dim<3> >*>&,
                                       const TableKernel<Dim<3> >&,
                                       const SmoothingScaleBase<Dim<3> >&,
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
#endif
}