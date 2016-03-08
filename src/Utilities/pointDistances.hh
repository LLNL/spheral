//---------------------------------Spheral++----------------------------------//
// pointDistances
//
// Assorted utilities for finding the distances from a point to things.
//
// Created by JMO, Mon Apr 11 21:19:13 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_pointDistances_hh__
#define __Spheral_pointDistances_hh__

#include "Geometry/Dimension.hh"
#include "safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Find the point on a line segment (a0,a1) closest to point (p).
//------------------------------------------------------------------------------
template<typename Vector>
inline
Vector
closestPointOnSegment(const Vector& p,
                      const Vector& a0,
                      const Vector& a1) {
  const Vector a01 = a1 - a0;
  const double a01mag = a01.magnitude();
  const Vector ahat = a01*safeInv(a01mag);
  const Vector a0p = p - a0;
  const double ptest = a0p.dot(ahat);
  if (ptest <= 0.0) return a0;
  if (ptest >= a01mag) return a1;
  return a0 + ptest*ahat;
}

//------------------------------------------------------------------------------
// Compute the distance between a point and a plane.
//------------------------------------------------------------------------------
template<typename Vector>
inline
double
pointPlaneDistance(const Vector& point,
                   const Vector& origin,
                   const Vector& unitNormal) {
   REQUIRE(fuzzyEqual(unitNormal.magnitude2(), 1.0, 1.0e-10));
   return std::abs((point - origin).dot(unitNormal));
}

//------------------------------------------------------------------------------
// Find the point on a plane closest to the given point.
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
closestPointOnPlane(const Dim<3>::Vector& p,
                    const Dim<3>::Vector& origin,
                    const Dim<3>::Vector& unitNormal) {
  REQUIRE(fuzzyEqual(unitNormal.magnitude2(), 1.0, 1.0e-10));
  return p + (origin - p).dot(unitNormal)*unitNormal;
}

}

#endif
