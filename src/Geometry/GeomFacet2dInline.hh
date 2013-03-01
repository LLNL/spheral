#include "GeomVector.hh"
#include "Infrastructure/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors, destructor.
//------------------------------------------------------------------------------
// We really don't want the default constructor, but it's required to have 
// std::vectors of these.
inline
GeomFacet2d::
GeomFacet2d():
  mVerticesPtr(0),
  mPoint1(0),
  mPoint2(0),
  mNormal(1.0, 0.0) {
  VERIFY(false);
}

inline
GeomFacet2d::
GeomFacet2d(const std::vector<GeomFacet2d::Vector>& vertices,
            const unsigned point1,
            const unsigned point2):
  mVerticesPtr(&vertices),
  mPoint1(point1),
  mPoint2(point2),
  mNormal(vertices[point2].y() - vertices[point1].y(),
          vertices[point1].x() - vertices[point2].x()) {
  REQUIRE((this->point2() - this->point1()).magnitude2() > 0.0);
  REQUIRE(fuzzyEqual((this->point2() - this->point1()).dot(mNormal), 0.0, 1.0e-10));
  REQUIRE((this->point2() - this->point1()).cross(mNormal).z() <= 0.0);
}

inline
GeomFacet2d::
GeomFacet2d(const GeomFacet2d& rhs):
  mVerticesPtr(rhs.mVerticesPtr),
  mPoint1(rhs.mPoint1),
  mPoint2(rhs.mPoint2),
  mNormal(rhs.mNormal) {
}

inline
GeomFacet2d&
GeomFacet2d::
operator=(const GeomFacet2d& rhs) {
  if (this != &rhs) {
    mVerticesPtr = rhs.mVerticesPtr;
    mPoint1 = rhs.mPoint1;
    mPoint2 = rhs.mPoint2;
    mNormal = rhs.mNormal;
  }
  return *this;
}

inline
GeomFacet2d::
~GeomFacet2d() {
}

//------------------------------------------------------------------------------
// Is the given point above, below, or colinear with the facet?
// Returns 1, -1, 0 respectively.
//------------------------------------------------------------------------------
inline
int
GeomFacet2d::
compare(const GeomFacet2d::Vector& point,
        const double tol) const {
  const double test = mNormal.dot(point - (*mVerticesPtr)[mPoint1]);
  return (fuzzyEqual(test, 0.0, tol) ?  0 :
          test > 0.0                 ?  1 :
                                       -1);
}

//------------------------------------------------------------------------------
// Access the points and normal.
//------------------------------------------------------------------------------
inline
const GeomFacet2d::Vector&
GeomFacet2d::
point1() const {
  REQUIRE(mVerticesPtr != 0 and
          mPoint1 < mVerticesPtr->size());
  return (*mVerticesPtr)[mPoint1];
}

inline
const GeomFacet2d::Vector&
GeomFacet2d::
point2() const {
  REQUIRE(mVerticesPtr != 0 and
          mPoint2 < mVerticesPtr->size());
  return (*mVerticesPtr)[mPoint2];
}

inline
unsigned
GeomFacet2d::
ipoint1() const {
  return mPoint1;
}

inline
unsigned
GeomFacet2d::
ipoint2() const {
  return mPoint2;
}

inline
const GeomFacet2d::Vector&
GeomFacet2d::
normal() const {
  return mNormal;
}

//------------------------------------------------------------------------------
// Compute the position.
//------------------------------------------------------------------------------
inline
GeomFacet2d::Vector
GeomFacet2d::
position() const {
  return 0.5*(point1() + point2());
}

//------------------------------------------------------------------------------
// Compute the area (line length in this case).
//------------------------------------------------------------------------------
inline
double
GeomFacet2d::
area() const {
  return (point2() - point1()).magnitude();
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
inline
bool
GeomFacet2d::
operator==(const GeomFacet2d& rhs) const {
  return (*mVerticesPtr == *(rhs.mVerticesPtr) and
          mPoint1 == rhs.mPoint1 and
          mPoint2 == rhs.mPoint2);
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
inline
bool
GeomFacet2d::
operator!=(const GeomFacet2d& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
inline
std::ostream&
operator<<(std::ostream& os, const GeomFacet2d& facet) {
  os << "GeomFacet2d( ivertices : " << facet.ipoint1() << " " << facet.ipoint2() << "\n"
     << "              vertices : " << facet.point1() << " " << facet.point2() << "\n"
     << "                normal : " << facet.normal() 
     << "\n)";
  return os;
}

}
