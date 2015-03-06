//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH interpolation at each point.
//------------------------------------------------------------------------------
#ifndef __Spheral__interpolateCRKSPH__
#define __Spheral__interpolateCRKSPH__

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

  namespace CRKSPHSpace {

    template<typename Dimension, typename DataType>
    FieldSpace::FieldList<Dimension, DataType>
    interpolateCRKSPH(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
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
