//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "iterateIdealH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void iterateIdealH<Dim<1> >(DataBase<Dim<1> >&, 
                                       const vector<Boundary<Dim<1> >*>&, 
                                       const TableKernel<Dim<1> >&,
                                       const SmoothingScaleBase<Dim<1> >&,
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
  template void iterateIdealH<Dim<2> >(DataBase<Dim<2> >&, 
                                       const vector<Boundary<Dim<2> >*>&, 
                                       const TableKernel<Dim<2> >&,
                                       const SmoothingScaleBase<Dim<2> >&,
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
  template void iterateIdealH<Dim<3> >(DataBase<Dim<3> >&, 
                                       const vector<Boundary<Dim<3> >*>&, 
                                       const TableKernel<Dim<3> >&,
                                       const SmoothingScaleBase<Dim<3> >&,
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
}
