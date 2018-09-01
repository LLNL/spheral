//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSolidCRKSPHSumMassDensity__
#define __Spheral__computeSolidCRKSPHSumMassDensity__

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class Boundary;
class NodeCoupling;

template<typename Dimension>
void
computeSolidCRKSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                 const KernelSpace::TableKernel<Dimension>& W,
                                 const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                                 const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                                 const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                                 const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity0,
                                 const NodeCoupling& nodeCoupling,
                                 FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);

}

#endif
