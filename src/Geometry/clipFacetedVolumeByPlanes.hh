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
// Created by J. Michael Owen, Tue Nov 28 10:00:51 PST 2017
//----------------------------------------------------------------------------//
#ifndef __Spheral_clipFacetedVolumeByPlane__
#define __Spheral_clipFacetedVolumeByPlane__

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

namespace Spheral {

// 1D
void clipFacetedVolumeByPlanes(Dim<1>::FacetedVolume& poly, const std::vector<GeomPlane<Dim<1>>>& planes);

// 3D
void clipFacetedVolumeByPlanes(Dim<3>::FacetedVolume& poly, const std::vector<GeomPlane<Dim<3>>>& planes);

}

#endif

