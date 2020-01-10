//---------------------------------Spheral++----------------------------------//
// Geom3Vector -- Explicitly 3-D Geometric Vector Class.
//
// Created by JNJ, Wed Dec 21 13:13:00 PST 2005
//----------------------------------------------------------------------------//

#ifndef __Spheral__Geom3Vector_hh__
#define __Spheral__Geom3Vector_hh__

#include "Geometry/GeomVector.hh"

namespace Spheral {

// Three-dimensional vectors can exist in any dimension, since they may not 
// interact with or create forces in transverse dimensions.  Magnetic field 
// lines in a conducting fluid are a good example of three-dimensional 
// vectors embedded in a lower-dimensional space.  We need only create a 
// type that is distinct from GeomVector<3>, which we do below.
class Geom3Vector
{
  public:

  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;

  // Useful static member data.
  static const size_type nDimensions;
  static const size_type numElements;
  static const Geom3Vector zero;

  // Constructors.
  explicit Geom3Vector(const double x = 0.0,
                       const double y = 0.0,
                       const double z = 0.0);

  Geom3Vector(const Geom3Vector& vec);

  // Implicit conversion from vectors that happen to be 3-D.
  Geom3Vector(const GeomVector<3>& vec);

  // Destructor.
  ~Geom3Vector();

  // Assignment operators.
  Geom3Vector& operator=(const Geom3Vector& rhs);

  // Access elements by indicies.
  double operator()(size_type index) const;
  double& operator()(size_type index);

  // Access the individual elements by (x, y, z) notation.
  double x() const;
  double y() const;
  double z() const;
  void x(const double val);
  void y(const double val);
  void z(const double val);

  // Iterator access to the raw data.
  iterator begin() { return mGeomVector.begin(); }
  iterator end() { return mGeomVector.end(); }

  const_iterator begin() const { return mGeomVector.begin(); }
  const_iterator end() const { return mGeomVector.end(); }

  // Zero the vector.
  void Zero() { mGeomVector.Zero(); }

  // Mathematical operators.
  Geom3Vector operator-() const { return Geom3Vector(-mGeomVector); }

  Geom3Vector operator+(const Geom3Vector& vec) const { return Geom3Vector(mGeomVector + vec.mGeomVector); }
  Geom3Vector operator-(const Geom3Vector& vec) const { return Geom3Vector(mGeomVector - vec.mGeomVector); }
  Geom3Vector operator*(const double val) const { return Geom3Vector(mGeomVector * val); }
  Geom3Vector operator/(const double val) const { return Geom3Vector(mGeomVector / val); }

  Geom3Vector& operator+=(const Geom3Vector& vec) { mGeomVector += vec.mGeomVector; return *this; }
  Geom3Vector& operator-=(const Geom3Vector& vec) { mGeomVector -= vec.mGeomVector; return *this; }
  Geom3Vector& operator*=(const double val) { mGeomVector *= val; return *this; }
  Geom3Vector& operator/=(const double val) { mGeomVector /= val; return *this; }

  bool operator==(const Geom3Vector& vec) const { return mGeomVector == vec.mGeomVector; }
  bool operator!=(const Geom3Vector& vec) const { return mGeomVector != vec.mGeomVector; }
  bool operator<(const Geom3Vector& vec) const { return mGeomVector < vec.mGeomVector; }
  bool operator>(const Geom3Vector& vec) const { return mGeomVector > vec.mGeomVector; }
  bool operator<=(const Geom3Vector& vec) const { return mGeomVector <= vec.mGeomVector; }
  bool operator>=(const Geom3Vector& vec) const { return mGeomVector >= vec.mGeomVector; }

  double dot(const Geom3Vector& rhs) const { return mGeomVector.dot(rhs.mGeomVector); }
  Geom3Vector cross(const Geom3Vector& rhs) const { return Geom3Vector(mGeomVector.cross(rhs.mGeomVector)); }
  GeomTensor<3> dyad(const Geom3Vector& rhs) const;
  GeomSymmetricTensor<3> selfdyad() const;
  GeomTensor<3> operator*(const Geom3Vector& rhs) const;

  Geom3Vector unitVector() const { return Geom3Vector(mGeomVector.unitVector()); }

  double magnitude() const { return mGeomVector.magnitude(); }
  double magnitude2() const { return mGeomVector.magnitude2(); }
  double minElement() const { return mGeomVector.minElement(); }
  double maxElement() const { return mGeomVector.maxElement(); }
  double sumElements() const { return mGeomVector.sumElements(); }

private:
  GeomVector<3> mGeomVector;

};

std::istream& operator>>(std::istream& is, Geom3Vector& vec);
std::ostream& operator<<(std::ostream& os, const Geom3Vector& vec);

}

#ifndef __GCCXML__
#include "Geom3VectorInline.hh"
#endif

#else

namespace Spheral {
  class Geom3Vector;
}

#endif

