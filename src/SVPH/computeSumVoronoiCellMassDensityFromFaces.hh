//------------------------------------------------------------------------------
// Compute the Voronoi cell mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeVoronoiCellMassDensityFromFaces__
#define __Spheral__computeVoronoiCellMassDensityFromFaces__

namespace Spheral {

  // Forward declarations.
  namespace MeshSpace {
    template<typename Dimension> class Mesh;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }

  namespace SVPHSpace {

    template<typename Dimension>
    void
    computeSumVoronoiCellMassDensityFromFaces(const MeshSpace::Mesh<Dimension>& mesh,
                                              const KernelSpace::TableKernel<Dimension>& W,
                                              const DataBaseSpace::DataBase<Dimension>& dataBase,
                                              FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);

  }
}

#endif
