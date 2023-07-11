#include "GeomVector.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Is the given point above, below, or coplanar with the facet?
// Returns 1, -1, 0 respectively.
//------------------------------------------------------------------------------
inline
int
GeomFacet3d::
compare(const GeomFacet3d::Vector& point,
        const double tol) const {
  const double test = mNormal.dot(point - (*mVerticesPtr)[mPoints[0]]);
  return (fuzzyEqual(test, 0.0, tol) ?  0 :
          test > 0.0                 ?  1 :
                                       -1);
}

//------------------------------------------------------------------------------
// Access the points and normal.
//------------------------------------------------------------------------------
inline
const GeomFacet3d::Vector&
GeomFacet3d::
point(const unsigned index) const {
  REQUIRE(mVerticesPtr != 0 and
          index < mPoints.size() and
          mPoints[index] < mVerticesPtr->size());
  return (*mVerticesPtr)[mPoints[index]];
}

inline
const std::vector<unsigned>&
GeomFacet3d::
ipoints() const {
  return mPoints;
}

inline
const GeomFacet3d::Vector&
GeomFacet3d::
normal() const {
  return mNormal;
}

}
