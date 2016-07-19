//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computecomputeVoronoiVolume__
#define __Spheral__computecomputeVoronoiVolume__

#include <vector>

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPHCorrectionParams.hh"

namespace Spheral {
  namespace CRKSPHSpace {

#ifdef SPHERAL1D
    // 1D
    void
    computeVoronoiVolume(const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& position,
                         const FieldSpace::FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                         const NeighborSpace::ConnectivityMap<Dim<1> >& connectivityMap,
                         const Dim<1>::Scalar kernelExtent,
                         FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& vol);
#endif

#ifdef SPHERAL2D
    // 2D
    void
    computeVoronoiVolume(const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& position,
                         const FieldSpace::FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                         const NeighborSpace::ConnectivityMap<Dim<2> >& connectivityMap,
                         const Dim<2>::Scalar kernelExtent,
                         FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>& vol);
#endif

#ifdef SPHERAL3D
    // 3D
    void
    computeVoronoiVolume(const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& position,
                         const FieldSpace::FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                         const NeighborSpace::ConnectivityMap<Dim<3> >& connectivityMap,
                         const Dim<3>::Scalar kernelExtent,
                         FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& vol);
#endif

  }
}

#endif
