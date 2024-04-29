text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/iterateIdealH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void iterateIdealH<Dim< %(ndim)s > >(DataBase<Dim< %(ndim)s > >&, 
                                       SmoothingScaleBase<Dim< %(ndim)s > >&,
                                       const vector<Boundary<Dim< %(ndim)s > >*>&, 
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
}
"""
