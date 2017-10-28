#include "GeomVector.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors, destructor.
//------------------------------------------------------------------------------
// We really don't want the default constructor, but it's required to have 
// std::vectors of these.
inline
GeomFacet3d::
GeomFacet3d():
  mVerticesPtr(0),
  mPoints(),
  mNormal(1.0, 0.0, 0.0) {
  VERIFY(false);
}

inline
GeomFacet3d::
GeomFacet3d(const std::vector<GeomFacet3d::Vector>& vertices,
            const std::vector<unsigned>& ipoints,
            const GeomFacet3d::Vector& normal):
  mVerticesPtr(&vertices),
  mPoints(ipoints),
  mNormal(normal) {
  REQUIRE(mPoints.size() >= 3);
}

inline
GeomFacet3d::
GeomFacet3d(const GeomFacet3d& rhs):
  mVerticesPtr(rhs.mVerticesPtr),
  mPoints(rhs.mPoints),
  mNormal(rhs.mNormal) {
}

inline
GeomFacet3d&
GeomFacet3d::
operator=(const GeomFacet3d& rhs) {
  if (this != &rhs) {
    mVerticesPtr = rhs.mVerticesPtr;
    mPoints = rhs.mPoints;
    mNormal = rhs.mNormal;
  }
  return *this;
}

inline
GeomFacet3d::
~GeomFacet3d() {
}

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

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
inline
bool
GeomFacet3d::
operator==(const GeomFacet3d& rhs) const {
  return (*mVerticesPtr == *(rhs.mVerticesPtr) and
          mPoints == rhs.mPoints and
          mNormal == rhs.mNormal);
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
inline
bool
GeomFacet3d::
operator!=(const GeomFacet3d& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
inline
std::ostream&
operator<<(std::ostream& os, const GeomFacet3d& facet) {
  os << "GeomFacet3d( ivertices : ";
  const std::vector<unsigned>& ipoints = facet.ipoints();
  for (unsigned i = 0; i != ipoints.size(); ++i) os << ipoints[i] << " ";
  os << "\n              vertices : ";
  for (unsigned i = 0; i != ipoints.size(); ++i) os << facet.point(i) << " ";
  os << "\n                normal : " << facet.normal() 
     << "\n)";
  return os;
}

}
