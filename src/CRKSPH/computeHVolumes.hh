//---------------------------------Spheral++------------------------------------
// Compute a volume per point based on the local H tensor.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeHVolumes__
#define __Spheral__computeHVolumes__

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

  template<typename Dimension>
  void
  computeHVolumes(const typename Dimension::Scalar nPerh,
                  const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                  FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& volume);
}

#endif
