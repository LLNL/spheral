//---------------------------------Spheral++------------------------------------
// Compute the TaylorSPH corrections.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeTaylorSPHCorrections__
#define __Spheral__computeTaylorSPHCorrections__

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
void
computeTaylorSPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Scalar>& weight,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Tensor>& D);

}

#endif
