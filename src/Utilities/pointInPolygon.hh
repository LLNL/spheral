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
#ifndef __Spheral_pointInPolygon__
#define __Spheral_pointInPolygon__

#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Work with a closed polygon defined by it's vertices.
// We assume the caller has ordered the vertices (clockwise or 
// counter-clockwise), and should have already performed the box exclusion
// test for efficiency.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<2>::Vector& p,
                    const std::vector<Dim<2>::Vector>& vertices);

//------------------------------------------------------------------------------
// Test a polygon (2-D).
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<2>::Vector& p,
                    const Dim<2>::FacetedVolume& polygon,
                    const bool countBoundary = false,
                    const double tol = 1.0e-10);

//------------------------------------------------------------------------------
// A 3-D version -- the point and all vertices of the polygon must be coplanar.
// p : point we're testing
// vertices : a set of 3-D positions, some subset of which define the polygon 
//            we're testing.
//
// The second form allows the user to specify a subset of the vertices to 
// represent the polygon as ipoints.
//------------------------------------------------------------------------------
bool pointInPolygon(const Dim<3>::Vector& p,
                    const std::vector<Dim<3>::Vector>& vertices,
                    const Dim<3>::Vector& normal);

bool pointInPolygon(const Dim<3>::Vector& p,
                    const std::vector<Dim<3>::Vector>& vertices,
                    const std::vector<unsigned>& ipoints,
                    const Dim<3>::Vector& normal,
                    const bool countBoundary = false,
                    const double tol = 1.0e-10);

}

#endif
