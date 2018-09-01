//---------------------------------Spheral++------------------------------------
// Compute the hull mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeHullSumMassDensity__
#define __Spheral__computeHullSumMassDensity__

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class Boundary;
class NodeCoupling;

template<typename Dimension>
void
computeHullSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                          const KernelSpace::TableKernel<Dimension>& W,
                          const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                          const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                          const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                          const NodeCoupling& nodeCoupling,
                          FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);

}

#endif
