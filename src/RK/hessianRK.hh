//---------------------------------Spheral++------------------------------------
// Compute the RK second-derivative
//------------------------------------------------------------------------------
#ifndef __Spheral__hessianRK__
#define __Spheral__hessianRK__

#include "Geometry/MathTraits.hh"
#include "Utilities/NodeCoupling.hh"
#include "RK/RKCorrectionParams.hh"
#include "RK/ReproducingKernel.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::HessianType>
hessianRK(const FieldList<Dimension, DataType>& fieldList,
          const FieldList<Dimension, typename Dimension::Vector>& position,
          const FieldList<Dimension, typename Dimension::Scalar>& weight,
          const FieldList<Dimension, typename Dimension::SymTensor>& H,
          const ConnectivityMap<Dimension>& connectivityMap,
          const ReproducingKernel<Dimension>& WR,
          const FieldList<Dimension, RKCoefficients<Dimension>>& corrections,
          const NodeCoupling& nodeCoupling = NodeCoupling());

}

#endif
