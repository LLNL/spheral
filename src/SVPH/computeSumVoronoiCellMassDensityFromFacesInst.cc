//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "computeSumVoronoiCellMassDensityFromFaces.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SVPHSpace {
    template void computeSumVoronoiCellMassDensityFromFaces(const Mesh<Dim<1> >&,
                                                            const TableKernel<Dim<1> >&, 
                                                            const DataBase<Dim<1> >&,
                                                            FieldList<Dim<1>, Dim<1>::Scalar>&);
    template void computeSumVoronoiCellMassDensityFromFaces(const Mesh<Dim<2> >&,
                                                            const TableKernel<Dim<2> >&, 
                                                            const DataBase<Dim<2> >&,
                                                            FieldList<Dim<2>, Dim<2>::Scalar>&);
    template void computeSumVoronoiCellMassDensityFromFaces(const Mesh<Dim<3> >&,
                                                            const TableKernel<Dim<3> >&, 
                                                            const DataBase<Dim<3> >&,
                                                            FieldList<Dim<3>, Dim<3>::Scalar>&);
  }
}

