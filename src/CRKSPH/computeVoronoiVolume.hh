//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computecomputeVoronoiVolume__
#define __Spheral__computecomputeVoronoiVolume__

#include <vector>

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {
  namespace CRKSPHSpace {

#ifdef SPHERAL1D
    // 1D
    void
    computeVoronoiVolume(const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& position,
                         const FieldSpace::FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                         const FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& rho,
                         const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& gradRho,
                         const NeighborSpace::ConnectivityMap<Dim<1> >& connectivityMap,
                         const std::vector<Dim<1>::FacetedVolume>& facetedBoundaries,
                         const std::vector<std::vector<Dim<1>::FacetedVolume> >& holes,
                         const std::vector<BoundarySpace::Boundary<Dim<1>>*>& boundaries,
                         const FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                         const FieldSpace::FieldList<Dim<1>, int>& voidPoint,
                         FieldSpace::FieldList<Dim<1>, int>& surfacePoint,
                         FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                         FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& deltaMedian,
                         FieldSpace::FieldList<Dim<1>, std::vector<Dim<1>::Vector>>& etaVoidPoints,
                         FieldSpace::FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells);
#endif

#ifdef SPHERAL2D
    // 2D
    void
    computeVoronoiVolume(const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& position,
                         const FieldSpace::FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                         const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>& rho,
                         const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& gradRho,
                         const NeighborSpace::ConnectivityMap<Dim<2> >& connectivityMap,
                         const std::vector<Dim<2>::FacetedVolume>& facetedBoundaries,
                         const std::vector<std::vector<Dim<2>::FacetedVolume> >& holes,
                         const std::vector<BoundarySpace::Boundary<Dim<2>>*>& boundaries,
                         const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                         const FieldSpace::FieldList<Dim<2>, int>& voidPoint,
                         FieldSpace::FieldList<Dim<2>, int>& surfacePoint,
                         FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>& vol,
                         FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& deltaMedian,
                         FieldSpace::FieldList<Dim<2>, std::vector<Dim<2>::Vector>>& etaVoidPoints,
                         FieldSpace::FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells);
#endif

#ifdef SPHERAL3D
    // 3D
    void
    computeVoronoiVolume(const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& position,
                         const FieldSpace::FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                         const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& rho,
                         const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& gradRho,
                         const NeighborSpace::ConnectivityMap<Dim<3> >& connectivityMap,
                         const std::vector<Dim<3>::FacetedVolume>& facetedBoundaries,
                         const std::vector<std::vector<Dim<3>::FacetedVolume> >& holes,
                         const std::vector<BoundarySpace::Boundary<Dim<3>>*>& boundaries,
                         const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                         const FieldSpace::FieldList<Dim<3>, int>& voidPoint,
                         FieldSpace::FieldList<Dim<3>, int>& surfacePoint,
                         FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& vol,
                         FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& deltaMedian,
                         FieldSpace::FieldList<Dim<3>, std::vector<Dim<3>::Vector>>& etaVoidPoints,
                         FieldSpace::FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells);
#endif

  }
}

#endif
