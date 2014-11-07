//---------------------------------Spheral++------------------------------------
// Compute the CSPH interpolation at each point.
//------------------------------------------------------------------------------
#ifndef __Spheral__interpolateCSPH__
#define __Spheral__interpolateCSPH__

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

  namespace CSPHSpace {

    template<typename Dimension, typename DataType>
    FieldSpace::FieldList<Dimension, DataType>
    interpolateCSPH(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
                    const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                    const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                    const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                    const bool coupleNodeLists,
                    const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                    const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                    const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                    const KernelSpace::TableKernel<Dimension>& W);

  }
}

#endif
