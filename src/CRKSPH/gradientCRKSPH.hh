//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH gradient.
//------------------------------------------------------------------------------
#ifndef __Spheral__gradientCRKSPH__
#define __Spheral__gradientCRKSPH__

#include "Geometry/MathTraits.hh"
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
    FieldSpace::FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
    gradientCRKSPH(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
                   const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                   const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                   const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                   const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                   const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                   const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& C,
                   const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradA,
                   const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradB,
                   const FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
                   const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
		   const int correctionOrder,
                   const KernelSpace::TableKernel<Dimension>& W,
                   const NodeCoupling& nodeCoupling = NodeCoupling());

  }
}

#endif
