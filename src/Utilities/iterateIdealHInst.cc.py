text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "iterateIdealH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void iterateIdealH<Dim< %(ndim)s > >(DataBase<Dim< %(ndim)s > >&, 
                                       const vector<Boundary<Dim< %(ndim)s > >*>&, 
                                       const TableKernel<Dim< %(ndim)s > >&,
                                       const SmoothingScaleBase<Dim< %(ndim)s > >&,
                                       const int,
                                       const double,
                                       const double,
                                       const bool,
                                       const bool);
}
"""
