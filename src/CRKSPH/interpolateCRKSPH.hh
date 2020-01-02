//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH interpolation at each point.
//------------------------------------------------------------------------------
#ifndef __Spheral__interpolateCRKSPH__
#define __Spheral__interpolateCRKSPH__

#include "SPH/NodeCoupling.hh"
#include "RK/RKCorrectionParams.hh"

#include "boost/variant.hpp"
#include <vector>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
std::vector<boost::variant<FieldList<Dimension, typename Dimension::Scalar>,
                           FieldList<Dimension, typename Dimension::Vector>,
                           FieldList<Dimension, typename Dimension::Tensor>,
                           FieldList<Dimension, typename Dimension::SymTensor>,
                           FieldList<Dimension, typename Dimension::ThirdRankTensor>>>
interpolateCRKSPH(const std::vector<boost::variant<FieldList<Dimension, typename Dimension::Scalar>,
                                                   FieldList<Dimension, typename Dimension::Vector>,
                                                   FieldList<Dimension, typename Dimension::Tensor>,
                                                   FieldList<Dimension, typename Dimension::SymTensor>,
                                                   FieldList<Dimension, typename Dimension::ThirdRankTensor>>>& fieldLists,
                  const FieldList<Dimension, typename Dimension::Vector>& position,
                  const FieldList<Dimension, typename Dimension::Scalar>& weight,
                  const FieldList<Dimension, typename Dimension::SymTensor>& H,
                  const FieldList<Dimension, typename Dimension::Scalar>& A,
                  const FieldList<Dimension, typename Dimension::Vector>& B,
                  const FieldList<Dimension, typename Dimension::Tensor>& C,
                  const ConnectivityMap<Dimension>& connectivityMap,
                  const RKOrder correctionOrder,
                  const TableKernel<Dimension>& W,
                  const NodeCoupling& nodeCoupling = NodeCoupling());

}

#endif
