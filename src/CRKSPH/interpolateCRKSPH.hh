//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH interpolation at each point.
//------------------------------------------------------------------------------
#ifndef __Spheral__interpolateCRKSPH__
#define __Spheral__interpolateCRKSPH__

#include "SolidSPH/NodeCoupling.hh"

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
                      const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& C,
                      const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                      const int correctionOrder,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const NodeCoupling& nodeCoupling = NodeCoupling());

  }
}

#endif
