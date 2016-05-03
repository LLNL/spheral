//------------------------------------------------------------------------------
// Use SVPH to take the gradient of a FieldList.
//------------------------------------------------------------------------------
#ifndef __Spheral__gradientFieldListSVPH__
#define __Spheral__gradientFieldListSVPH__

#include <vector>
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Mesh/Mesh.hh"
#include "Geometry/MathTraits.hh"

namespace Spheral {
  namespace SVPHSpace {

    template<typename Dimension, typename DataType>
    FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
    gradientFieldListSVPH(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
                          const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                          const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                          const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                          const KernelSpace::TableKernel<Dimension>& W,
                          const MeshSpace::Mesh<Dimension>& mesh,
                          const bool firstOrderConsistent);
  }
}

#endif
