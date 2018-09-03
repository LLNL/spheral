//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on convex hulls.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeHullVolumes__
#define __Spheral__computeHullVolumes__

namespace Spheral {

// Forward declarations.
template<typename Dimension> class ConnectivityMap;
template<typename Dimension> class TableKernel;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
void
computeHullVolumes(const ConnectivityMap<Dimension>& connectivityMap,
                   const typename Dimension::Scalar kernelExtent,
                   const FieldList<Dimension, typename Dimension::Vector>& position,
                   const FieldList<Dimension, typename Dimension::SymTensor>& H,
                   FieldList<Dimension, typename Dimension::Scalar>& volume);

}

#endif
