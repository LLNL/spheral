//---------------------------------Spheral++----------------------------------//
// Clip a faceted volume (box, polygon, or polyhedron) by a plane in place.
//
// We use the convention that any portion of the faceted volume "below" the 
// plane is clipped, i.e., only the portion of the faceted volume "above" the 
// plane
//    plane.compare(point) >= 0
// is retained.
//
// The algorithms herein are roughly based on the approach outlined in 
// Powell, D., & Abel, T. (2015). An exact general remeshing scheme applied to 
// physically conservative voxelization. Journal of Computational Physics, 297, 340–356.
// though I think not exactly the same.
//
// This version is specialized for convex polygons/polyhedra.
//
// Created by J. Michael Owen, Wed Dec 13 15:28:09 PST 2017
//----------------------------------------------------------------------------//
#ifndef __Spheral_clipConvexFacetedVolumeByPlane__
#define __Spheral_clipConvexFacetedVolumeByPlane__

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

namespace Spheral {

// 1D
void clipConvexFacetedVolumeByPlanes(Dim<1>::FacetedVolume& poly, const std::vector<GeomPlane<Dim<1>>>& planes);

// 2D
void clipConvexFacetedVolumeByPlanes(Dim<2>::FacetedVolume& poly, const std::vector<GeomPlane<Dim<2>>>& planes);

// 3D
void clipConvexFacetedVolumeByPlanes(Dim<3>::FacetedVolume& poly, const std::vector<GeomPlane<Dim<3>>>& planes);

}

#endif

