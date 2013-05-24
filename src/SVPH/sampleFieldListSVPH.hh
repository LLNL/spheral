//------------------------------------------------------------------------------
// Use SVPH to sample a FieldList.
//------------------------------------------------------------------------------
#ifndef __Spheral__sampleFieldListSVPH__
#define __Spheral__sampleFieldListSVPH__

#include "Geometry/MathTraits.hh"

namespace Spheral {

  // Forward declarations.
  namespace NeighborSpace {
    template<typename Dimension> class ConnectivityMap;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace MeshSpace {
    template<typename Dimension> class Mesh;
  }

  namespace SVPHSpace {

    template<typename Dimension, typename DataType>
    FieldList<Dimension, typename MathTraits<DataType>::GradientType>
    sampleFieldListSVPH(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                        const KernelSpace::TableKernel<Dimension>& W,
                        const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                        const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                        const MeshSpace::Mesh<Dimension>& mesh);
  }
}

#endif
