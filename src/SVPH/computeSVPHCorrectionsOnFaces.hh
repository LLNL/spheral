//------------------------------------------------------------------------------
// Modified ideas from CSPH using the SVPH normalized kernel estimates.
// This version computes the corrections on the faces of a mesh, rather than on
// the points themselves.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSVPHCorrectionsOnFaces__
#define __Spheral__computeSVPHCorrectionsOnFaces__

#include <vector>
#include "Geometry/Dimension.hh"

namespace Spheral {

  // Forward declarations.
  namespace MeshSpace {
    template<typename Dimension> class Mesh;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }

  namespace SVPHSpace {
    template<typename Dimension, typename BoundaryIterator>
    void
    computeSVPHCorrectionsOnFaces(const MeshSpace::Mesh<Dimension>& mesh,
                                  const KernelSpace::TableKernel<Dimension>& W,
                                  const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& volume,
                                  const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                                  const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                                  const BoundaryIterator& boundaryBegin,
                                  const BoundaryIterator& boundaryEnd,
                                  std::vector<typename Dimension::Scalar>& A,
                                  std::vector<typename Dimension::Vector>& B);
  }
}

#endif
