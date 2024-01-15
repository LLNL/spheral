//---------------------------------Spheral++----------------------------------//
// GeomPlane -- Geometric Plane Class.
//
// Created by JMO, Thu Feb 24 17:26:32 PST 2000
//----------------------------------------------------------------------------//

#ifndef __Spheral_GeomPlane__
#define __Spheral_GeomPlane__

#include <string>
#include <vector>
#include <iostream>

namespace Spheral {

template<typename Dimension>
class GeomPlane {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  GeomPlane();
  GeomPlane(const GeomPlane& rhs);
  GeomPlane(const Vector& point, const Vector& normal);

  // This special constructor finds a plane fitted to a collection of points.
  explicit GeomPlane(const std::vector<Vector>& points);

  // Destructor.
  ~GeomPlane();

  // Assignment.
  GeomPlane& operator=(const GeomPlane& rhs);

  // Access the point and normal which define the plane.
  const Vector& point() const;
  void point(const Vector& point);

  const Vector& normal() const;
  void normal(const Vector& normal);

  // Negative operator (reverses sign of normal).
  GeomPlane operator-() const;

  // Calculate the signed distance between a given point and the plane.
  double signedDistance(const Vector& point) const;

  // Calculate the minimum distance between a given point and the plane.
  double minimumDistance(const Vector& point) const;

  // Closest point on the plane to a point.
  Vector closestPointOnPlane(const Vector& p) const;

  // Various tests which can be applied between planes.
  bool parallel(const GeomPlane& rhs) const;
  bool operator==(const GeomPlane& rhs) const;
  bool operator!=(const GeomPlane& rhs) const;
  bool operator<(const GeomPlane& rhs) const;

  // Some tests which can applied between a plane and a point.
  // These are only meant as tests for whether a point lies "above" or "below" the
  // plane, as represented by the direction of the normal.
  bool operator==(const Vector& point) const;
  bool operator!=(const Vector& point) const;
  bool operator>(const Vector& point) const;
  bool operator<(const Vector& point) const;
  bool operator>=(const Vector& point) const;
  bool operator<=(const Vector& point) const;

  // Test if the plane is below, equal to or above the point (-1, 0, 1).
  int compare(const Vector& point) const;

  // The usual validity test.
  bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  Vector mPoint;
  Vector mNormal;
};

// Ostream operators.
template<typename Dimension> std::istream& operator>>(std::istream& is, GeomPlane<Dimension>& vec);
template<typename Dimension> std::ostream& operator<<(std::ostream& os, const GeomPlane<Dimension>& vec);

}

#include "GeomPlaneInline.hh"
#include "Dimension.hh"

// Declare specialized methods.
namespace Spheral {
template<> GeomPlane<Dim<1> >::GeomPlane(const std::vector<Dim<1>::Vector>& points);
template<> GeomPlane<Dim<2> >::GeomPlane(const std::vector<Dim<2>::Vector>& points);
template<> GeomPlane<Dim<3> >::GeomPlane(const std::vector<Dim<3>::Vector>& points);
}

#endif
