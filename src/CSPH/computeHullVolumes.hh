//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on convex hulls.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeHullVolumes__
#define __Spheral__computeHullVolumes__

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
  computeHullVolumes(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                     const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                     FieldSpace::FieldList<Dimension, typename Dimension::FacetedVolume>& polvol,
                     FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& volume);
}

#endif
