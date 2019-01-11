//------------------------------------------------------------------------------
// Compute the Voronoi cell mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeVoronoiCellMassDensityFromFaces__
#define __Spheral__computeVoronoiCellMassDensityFromFaces__

namespace Spheral {

// Forward declarations.
template<typename Dimension> class Mesh;
template<typename Dimension> class DataBase;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
void
computeSumVoronoiCellMassDensityFromFaces(const Mesh<Dimension>& mesh,
                                          const TableKernel<Dimension>& W,
                                          const DataBase<Dimension>& dataBase,
                                          FieldList<Dimension, typename Dimension::Scalar>& massDensity);

}

#endif
