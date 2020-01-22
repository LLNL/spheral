//---------------------------------Spheral++------------------------------------
// Compute a volume per point based on the local H tensor.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeHVolumes__
#define __Spheral__computeHVolumes__

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
void
computeHVolumes(const typename Dimension::Scalar nPerh,
                const FieldList<Dimension, typename Dimension::SymTensor>& H,
                FieldList<Dimension, typename Dimension::Scalar>& volume);

}

#endif
