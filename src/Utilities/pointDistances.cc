//---------------------------------Spheral++----------------------------------//
// pointDistances
//
// Assorted utilities for finding the distances from a point to things.
//
// Created by JMO, Mon Apr 11 21:19:13 PDT 2011
//----------------------------------------------------------------------------//
#include "pointDistances.hh"
#include "spheralWildMagicConverters.hh"

#include "Wm5Triangle3.h"
#include "Wm5DistPoint3Triangle3.h"

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the distance between a point (p) and a triangle (a, b, c).
//------------------------------------------------------------------------------
double
pointTriangleDistance(const Dim<3>::Vector& p,
                      const Dim<3>::Vector& a,
                      const Dim<3>::Vector& b,
                      const Dim<3>::Vector& c) {
  return Wm5::DistPoint3Triangle3<double>(convertVectorToWMVector<Dim<3> >(p),
                                          Wm5::Triangle3<double>(convertVectorToWMVector<Dim<3> >(a),
                                                                 convertVectorToWMVector<Dim<3> >(b),
                                                                 convertVectorToWMVector<Dim<3> >(c))).Get();
}

//------------------------------------------------------------------------------
// For a point (p) and a triangle (a, b, c), compute the coordinates of the 
// point on the triangle closest to p.
//------------------------------------------------------------------------------
Dim<3>::Vector
closestPointOnTriangle(const Dim<3>::Vector& p,
                       const Dim<3>::Vector& a,
                       const Dim<3>::Vector& b,
                       const Dim<3>::Vector& c) {
  const Wm5::Vector3<double> p_wm5 = convertVectorToWMVector<Dim<3> >(p);
  const Wm5::Triangle3<double> tri_wm5(convertVectorToWMVector<Dim<3> >(a),
                                       convertVectorToWMVector<Dim<3> >(b),
                                       convertVectorToWMVector<Dim<3> >(c));
  const Wm5::DistPoint3Triangle3<double> dist(p_wm5, tri_wm5);
  return convertWMVectorToVector<Dim<3> >(dist.GetTriangleBary(0)*tri_wm5.V[0] +
                                          dist.GetTriangleBary(1)*tri_wm5.V[1] +
                                          dist.GetTriangleBary(2)*tri_wm5.V[2]);
}

}
