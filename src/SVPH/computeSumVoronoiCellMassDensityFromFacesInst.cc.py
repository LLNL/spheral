text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/computeSumVoronoiCellMassDensityFromFaces.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template void computeSumVoronoiCellMassDensityFromFaces(const Mesh<Dim< %(ndim)s > >&,
                                                          const TableKernel<Dim< %(ndim)s > >&, 
                                                          const DataBase<Dim< %(ndim)s > >&,
                                                          FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}

"""
