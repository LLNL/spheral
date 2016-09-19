//------------------------------------------------------------------------------
// Compute the moments necessary for CRKSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHMoments__
#define __Spheral__computeCRKSPHMoments__

#include "SolidSPH/NodeCoupling.hh"
#include "CRKSPHCorrectionParams.hh"

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
    template<typename Dimension>
    void
    computeCRKSPHMoments(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                         const KernelSpace::TableKernel<Dimension>& W,
                         const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight,
                         const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                         const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                         const CRKOrder correctionOrder,
                         const NodeCoupling& nodeCoupling,
                         FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& m0,
                         FieldSpace::FieldList<Dimension, typename Dimension::Vector>& m1,
                         FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& m2,
                         FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& m3,
                         FieldSpace::FieldList<Dimension, typename Dimension::FourthRankTensor>& m4,
                         FieldSpace::FieldList<Dimension, typename Dimension::Vector>& gradm0,
                         FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                         FieldSpace::FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                         FieldSpace::FieldList<Dimension, typename Dimension::FourthRankTensor>& gradm3,
                         FieldSpace::FieldList<Dimension, typename Dimension::FifthRankTensor>& gradm4);
  }
}

#endif
