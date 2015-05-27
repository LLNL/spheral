//------------------------------------------------------------------------------
// pointInPolygon
//
// Test if a given point is in a polygon or not.  This is based on code I found
// at
//
// http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//
// which I have modified here somewhat to use Spheral structures and such.
//------------------------------------------------------------------------------
#ifndef __Spheral_pointInPolyhedron__
#define __Spheral_pointInPolyhedron__

#include "Geometry/Dimension.hh"
#include "testBoxIntersection.hh"
#include "pointOnPolygon.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Determine if a point is inside a polyhedron.
//------------------------------------------------------------------------------
bool pointInPolyhedron(const Dim<3>::Vector& p,
                       const Dim<3>::FacetedVolume& polyhedron,
                       const bool countBoundary = false,
                       const double tol = 1.0e-10);

}

#endif
