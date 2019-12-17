//------------------------------------------------------------------------------
// Compute the moments necessary for CRKSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHMoments__
#define __Spheral__computeCRKSPHMoments__

#include "SPH/NodeCoupling.hh"
#include "RK/RKCorrectionParams.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
void
computeCRKSPHMoments(const ConnectivityMap<Dimension>& connectivityMap,
                     const TableKernel<Dimension>& W,
                     const FieldList<Dimension, typename Dimension::Scalar>& weight,
                     const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     const RKOrder correctionOrder,
                     const NodeCoupling& nodeCoupling,
                     FieldList<Dimension, typename Dimension::Scalar>& m0,
                     FieldList<Dimension, typename Dimension::Vector>& m1,
                     FieldList<Dimension, typename Dimension::SymTensor>& m2,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& m3,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& m4,
                     FieldList<Dimension, typename Dimension::Vector>& gradm0,
                     FieldList<Dimension, typename Dimension::Tensor>& gradm1,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradm2,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& gradm3,
                     FieldList<Dimension, typename Dimension::FifthRankTensor>& gradm4);

}

#endif
