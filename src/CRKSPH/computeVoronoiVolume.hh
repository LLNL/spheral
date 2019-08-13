//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computecomputeVoronoiVolume__
#define __Spheral__computecomputeVoronoiVolume__

#include "Geometry/CellFaceFlag.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"

#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
#ifdef SPHERAL1D
// 1D
void
computeVoronoiVolume(const FieldList<Dim<1>, Dim<1>::Vector>& position,
                     const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                     const ConnectivityMap<Dim<1> >& connectivityMap,
                     const FieldList<Dim<1>, Dim<1>::SymTensor>& damage,
                     const std::vector<Dim<1>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<1>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<1>>*>& boundaries,
                     const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                     FieldList<Dim<1>, int>& surfacePoint,
                     FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                     FieldList<Dim<1>, Dim<1>::Vector>& deltaMedian,
                     FieldList<Dim<1>, std::vector<Dim<1>::Vector>>& etaVoidPoints,
                     FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells,
                     FieldList<Dim<1>, std::vector<CellFaceFlag>>& cellFaceFlags);
#endif

#ifdef SPHERAL2D
// 2D
void
computeVoronoiVolume(const FieldList<Dim<2>, Dim<2>::Vector>& position,
                     const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                     const ConnectivityMap<Dim<2> >& connectivityMap,
                     const FieldList<Dim<2>, Dim<2>::SymTensor>& damage,
                     const std::vector<Dim<2>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<2>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<2>>*>& boundaries,
                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                     FieldList<Dim<2>, int>& surfacePoint,
                     FieldList<Dim<2>, Dim<2>::Scalar>& vol,
                     FieldList<Dim<2>, Dim<2>::Vector>& deltaMedian,
                     FieldList<Dim<2>, std::vector<Dim<2>::Vector>>& etaVoidPoints,
                     FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells,
                     FieldList<Dim<2>, std::vector<CellFaceFlag>>& cellFaceFlags);
#endif

#ifdef SPHERAL3D
// 3D
void
computeVoronoiVolume(const FieldList<Dim<3>, Dim<3>::Vector>& position,
                     const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                     const ConnectivityMap<Dim<3> >& connectivityMap,
                     const FieldList<Dim<3>, Dim<3>::SymTensor>& damage,
                     const std::vector<Dim<3>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<3>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<3>>*>& boundaries,
                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                     FieldList<Dim<3>, int>& surfacePoint,
                     FieldList<Dim<3>, Dim<3>::Scalar>& vol,
                     FieldList<Dim<3>, Dim<3>::Vector>& deltaMedian,
                     FieldList<Dim<3>, std::vector<Dim<3>::Vector>>& etaVoidPoints,
                     FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells,
                     FieldList<Dim<3>, std::vector<CellFaceFlag>>& cellFaceFlags);
#endif

}

#endif
