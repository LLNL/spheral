//---------------------------------Spheral++----------------------------------//
// GeomPlane -- Geometric Plane Class.
//
// Created by JMO, Thu Feb 24 17:26:32 PST 2000
//----------------------------------------------------------------------------//

#ifndef GeomPlane_HH
#define GeomPlane_HH

#include <string>
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

  // Calculate the minimum distance between a given point and the plane.
  double minimumDistance(const Vector& point) const;

  // Various tests which can be applied between planes.
  bool parallel(const GeomPlane& rhs) const;
  bool operator==(const GeomPlane& rhs) const;
  bool operator!=(const GeomPlane& rhs) const;

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

template<typename Dimension> std::istream& operator>>(std::istream& is, GeomPlane<Dimension>& vec);
template<typename Dimension> std::ostream& operator<<(std::ostream& os, const GeomPlane<Dimension>& vec);

}

#ifndef __GCCXML__
#include "GeomPlaneInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class GeomPlane;
}

#endif
