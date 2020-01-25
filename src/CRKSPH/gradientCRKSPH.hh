//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH gradient.
//------------------------------------------------------------------------------
#ifndef __Spheral__gradientCRKSPH__
#define __Spheral__gradientCRKSPH__

#include "Geometry/MathTraits.hh"
#include "SPH/NodeCoupling.hh"
#include "RK/RKCorrectionParams.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientCRKSPH(const FieldList<Dimension, DataType>& fieldList,
               const FieldList<Dimension, typename Dimension::Vector>& position,
               const FieldList<Dimension, typename Dimension::Scalar>& weight,
               const FieldList<Dimension, typename Dimension::SymTensor>& H,
               const FieldList<Dimension, typename Dimension::Scalar>& A,
               const FieldList<Dimension, typename Dimension::Vector>& B,
               const FieldList<Dimension, typename Dimension::Tensor>& C,
               const FieldList<Dimension, typename Dimension::Vector>& gradA,
               const FieldList<Dimension, typename Dimension::Tensor>& gradB,
               const FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
               const ConnectivityMap<Dimension>& connectivityMap,
               const RKOrder correctionOrder,
               const TableKernel<Dimension>& W,
               const NodeCoupling& nodeCoupling = NodeCoupling());

}

#endif
