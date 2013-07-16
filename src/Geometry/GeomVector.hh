//---------------------------------Spheral++----------------------------------//
// GeomVector -- Geometric Vector Class.
//
// Created by JMO, Thu Apr 15 17:39:47 PDT 1999
// Modified by:
//   July 31, 99: JMO, removing MTL based version due to problems compiling
//                MTL code with gcc2.95.
//   2004-08-19:  JMO, switch to using boost::ublas, and our vector class is
//                just an interface.
//   2004-08-23:  JMO, ublas is still too slow, so going to primitive C 
//                internal data types in accordance with suggestions from
//                Brian White
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomVector_hh__
#define __Spheral_GeomVector_hh__

#include <iostream>

namespace Spheral {

template<int nDim> class GeomTensor;
template<int nDim> class GeomSymmetricTensor;

template<int nDim, bool ownMemory = true>
class GeomVector {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;

  // Useful static member data.
  static const int nDimensions;
  static const GeomVector zero;
  static const GeomVector one;

  // Constructors.
  GeomVector(const double x = 0.0,
             const double y = 0.0,
             const double z = 0.0);
  template<bool otherMemory> GeomVector(const GeomVector<nDim, otherMemory>& vec);

  // Destructor.
  ~GeomVector();

  // Assignment.
  template<bool otherMemory> GeomVector& operator=(const GeomVector<nDim, otherMemory>& vec);
  GeomVector& operator=(const double val);

  // Allow the elements by indicies.
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
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // Zero the vector.
  void Zero();

  // Mathematical operators.
  GeomVector<nDim, true> operator-() const;

  template<bool otherMemory> GeomVector<nDim, true> operator+(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> GeomVector<nDim, true> operator-(const GeomVector<nDim, otherMemory>& vec) const;
  GeomVector<nDim, true> operator+(const double val) const;
  GeomVector<nDim, true> operator-(const double val) const;
  GeomVector<nDim, true> operator*(const double val) const;
  GeomVector<nDim, true> operator/(const double val) const;

  template<bool otherMemory> GeomVector& operator+=(const GeomVector<nDim, otherMemory>& vec);
  template<bool otherMemory> GeomVector& operator-=(const GeomVector<nDim, otherMemory>& vec);
  GeomVector& operator+=(const double val);
  GeomVector& operator-=(const double val);
  GeomVector& operator*=(const double val);
  GeomVector& operator/=(const double val);

  template<bool otherMemory> int compare(const GeomVector<nDim, otherMemory>& vec) const;
  int compare(const double val) const;

  template<bool otherMemory> bool operator==(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator!=(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator< (const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator> (const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator<=(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> bool operator>=(const GeomVector<nDim, otherMemory>& vec) const;

  bool operator==(const double val) const;
  bool operator!=(const double val) const;
  bool operator< (const double val) const;
  bool operator> (const double val) const;
  bool operator<=(const double val) const;
  bool operator>=(const double val) const;

  template<bool otherMemory> double dot(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> GeomVector<3, true> cross(const GeomVector<nDim, otherMemory>& vec) const;
  template<bool otherMemory> GeomTensor<nDim> dyad(const GeomVector<nDim, otherMemory>& rhs) const;
  GeomSymmetricTensor<nDim> selfdyad() const;
  template<bool otherMemory> GeomTensor<nDim> operator*(const GeomVector<nDim, otherMemory>& vec) const;

  GeomVector<nDim, true> unitVector() const;

  double magnitude() const;
  double magnitude2() const;
  double minElement() const;
  double maxElement() const;
  double maxAbsElement() const;
  double sumElements() const;

private:
  double* mData;

  //  friend class GeomVector<nDim, !ownMemory>;
};

// Declare explicit specializations.
template<> GeomVector<1, true>::GeomVector(const double, const double, const double);
template<> GeomVector<2, true>::GeomVector(const double, const double, const double);
template<> GeomVector<3, true>::GeomVector(const double, const double, const double);

template<int nDim>
template<bool otherMemory>
GeomVector<nDim, true>::GeomVector(const GeomVector<nDim, otherMemory>& vec);

template<int nDim> GeomVector<nDim, true>::~GeomVector();
template<int nDim> GeomVector<nDim, false>::~GeomVector();

template<bool ownMemory> double GeomVector<1, ownMemory>::y() const;
template<bool ownMemory> double GeomVector<1, ownMemory>::z() const;
template<bool ownMemory> double GeomVector<2, ownMemory>::z() const;

template<bool ownMemory> void GeomVector<1, ownMemory>::y(const double val);
template<bool ownMemory> void GeomVector<1, ownMemory>::z(const double val);
template<bool ownMemory> void GeomVector<2, ownMemory>::z(const double val);

template<bool ownMemory> GeomVector<1, ownMemory> GeomVector<1, ownMemory>::operator-() const;
template<bool ownMemory> GeomVector<2, ownMemory> GeomVector<2, ownMemory>::operator-() const;
template<bool ownMemory> GeomVector<3, ownMemory> GeomVector<3, ownMemory>::operator-() const;

template<bool ownMemory> GeomVector<1, ownMemory>& GeomVector<1, ownMemory>::operator+=(const GeomVector<1, ownMemory>& vec);
template<bool ownMemory> GeomVector<2, ownMemory>& GeomVector<2, ownMemory>::operator+=(const GeomVector<2, ownMemory>& vec);
template<bool ownMemory> GeomVector<3, ownMemory>& GeomVector<3, ownMemory>::operator+=(const GeomVector<3, ownMemory>& vec);

template<bool ownMemory> GeomVector<1, ownMemory>& GeomVector<1, ownMemory>::operator-=(const GeomVector<1, ownMemory>& vec);
template<bool ownMemory> GeomVector<2, ownMemory>& GeomVector<2, ownMemory>::operator-=(const GeomVector<2, ownMemory>& vec);
template<bool ownMemory> GeomVector<3, ownMemory>& GeomVector<3, ownMemory>::operator-=(const GeomVector<3, ownMemory>& vec);

template<bool ownMemory> GeomVector<1, ownMemory>& GeomVector<1, ownMemory>::operator+=(const double val);
template<bool ownMemory> GeomVector<2, ownMemory>& GeomVector<2, ownMemory>::operator+=(const double val);
template<bool ownMemory> GeomVector<3, ownMemory>& GeomVector<3, ownMemory>::operator+=(const double val);

template<bool ownMemory> GeomVector<1, ownMemory>& GeomVector<1, ownMemory>::operator-=(const double val);
template<bool ownMemory> GeomVector<2, ownMemory>& GeomVector<2, ownMemory>::operator-=(const double val);
template<bool ownMemory> GeomVector<3, ownMemory>& GeomVector<3, ownMemory>::operator-=(const double val);

template<bool ownMemory> GeomVector<1, ownMemory>& GeomVector<1, ownMemory>::operator*=(const double val);
template<bool ownMemory> GeomVector<2, ownMemory>& GeomVector<2, ownMemory>::operator*=(const double val);
template<bool ownMemory> GeomVector<3, ownMemory>& GeomVector<3, ownMemory>::operator*=(const double val);

template<bool ownMemory> GeomVector<1, ownMemory>& GeomVector<1, ownMemory>::operator/=(const double val);
template<bool ownMemory> GeomVector<2, ownMemory>& GeomVector<2, ownMemory>::operator/=(const double val);
template<bool ownMemory> GeomVector<3, ownMemory>& GeomVector<3, ownMemory>::operator/=(const double val);

template<bool ownMemory> template<bool otherMemory> int GeomVector<1, ownMemory>::compare(const GeomVector<1, otherMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> int GeomVector<2, ownMemory>::compare(const GeomVector<2, otherMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> int GeomVector<3, ownMemory>::compare(const GeomVector<3, otherMemory>& vec) const;

template<bool ownMemory> int GeomVector<1, ownMemory>::compare(const double val) const;
template<bool ownMemory> int GeomVector<2, ownMemory>::compare(const double val) const;
template<bool ownMemory> int GeomVector<3, ownMemory>::compare(const double val) const;

template<bool ownMemory> template<bool otherMemory> bool GeomVector<1, ownMemory>::operator==(const GeomVector<1, ownMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> bool GeomVector<2, ownMemory>::operator==(const GeomVector<2, ownMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> bool GeomVector<3, ownMemory>::operator==(const GeomVector<3, ownMemory>& vec) const;

template<bool ownMemory> bool GeomVector<1, ownMemory>::operator==(const double val) const;
template<bool ownMemory> bool GeomVector<2, ownMemory>::operator==(const double val) const;
template<bool ownMemory> bool GeomVector<3, ownMemory>::operator==(const double val) const;

template<bool ownMemory> template<bool otherMemory> double GeomVector<1, ownMemory>::dot(const GeomVector<1, otherMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> double GeomVector<2, ownMemory>::dot(const GeomVector<2, otherMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> double GeomVector<3, ownMemory>::dot(const GeomVector<3, otherMemory>& vec) const;

template<bool ownMemory> template<bool otherMemory> GeomVector<3, true> GeomVector<1, ownMemory>::cross(const GeomVector<1, otherMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> GeomVector<3, true> GeomVector<2, ownMemory>::cross(const GeomVector<2, otherMemory>& vec) const;
template<bool ownMemory> template<bool otherMemory> GeomVector<3, true> GeomVector<3, ownMemory>::cross(const GeomVector<3, otherMemory>& vec) const;

template<bool ownMemory> template<bool otherMemory> GeomTensor<1, ownMemory> GeomVector<1, ownMemory>::dyad(const GeomVector<1, ownMemory>& rhs) const;
template<bool ownMemory> template<bool otherMemory> GeomTensor<2, ownMemory> GeomVector<2, ownMemory>::dyad(const GeomVector<2, ownMemory>& rhs) const;
template<bool ownMemory> template<bool otherMemory> GeomTensor<3, ownMemory> GeomVector<3, ownMemory>::dyad(const GeomVector<3, ownMemory>& rhs) const;

template<bool ownMemory> GeomSymmetricTensor<1, ownMemory> GeomVector<1, ownMemory>::selfdyad() const;
template<bool ownMemory> GeomSymmetricTensor<2, ownMemory> GeomVector<2, ownMemory>::selfdyad() const;
template<bool ownMemory> GeomSymmetricTensor<3, ownMemory> GeomVector<3, ownMemory>::selfdyad() const;

template<bool ownMemory> double GeomVector<1, ownMemory>::magnitude() const;
template<bool ownMemory> double GeomVector<2, ownMemory>::magnitude() const;
template<bool ownMemory> double GeomVector<3, ownMemory>::magnitude() const;

template<bool ownMemory> double GeomVector<1, ownMemory>::magnitude2() const;
template<bool ownMemory> double GeomVector<2, ownMemory>::magnitude2() const;
template<bool ownMemory> double GeomVector<3, ownMemory>::magnitude2() const;

template<bool ownMemory> double GeomVector<1, ownMemory>::minElement() const;
template<bool ownMemory> double GeomVector<2, ownMemory>::minElement() const;
template<bool ownMemory> double GeomVector<3, ownMemory>::minElement() const;

template<bool ownMemory> double GeomVector<1, ownMemory>::maxElement() const;
template<bool ownMemory> double GeomVector<2, ownMemory>::maxElement() const;
template<bool ownMemory> double GeomVector<3, ownMemory>::maxElement() const;

template<bool ownMemory> double GeomVector<1, ownMemory>::maxAbsElement() const;
template<bool ownMemory> double GeomVector<2, ownMemory>::maxAbsElement() const;
template<bool ownMemory> double GeomVector<3, ownMemory>::maxAbsElement() const;

template<bool ownMemory> double GeomVector<1, ownMemory>::sumElements() const;
template<bool ownMemory> double GeomVector<2, ownMemory>::sumElements() const;
template<bool ownMemory> double GeomVector<3, ownMemory>::sumElements() const;

// Forward declare the global functions.
template<int nDim, bool ownMemory> GeomVector<nDim> operator+(const double val, const GeomVector<nDim, ownMemory>& vec);
template<int nDim, bool ownMemory> GeomVector<nDim> operator-(const double val, const GeomVector<nDim, ownMemory>& vec);
template<int nDim, bool ownMemory> GeomVector<nDim> operator*(const double val, const GeomVector<nDim, ownMemory>& vec);

template<int nDim> GeomVector<nDim> elementWiseMin(const GeomVector<nDim>& lhs,
                                                   const GeomVector<nDim>& rhs);
template<int nDim> GeomVector<nDim> elementWiseMax(const GeomVector<nDim>& lhs,
                                                   const GeomVector<nDim>& rhs);

template<bool ownMemory, bool otherMemory> GeomVector<1> elementWiseMin(const GeomVector<1, ownMemory>& lhs,
                                                                        const GeomVector<1, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<2> elementWiseMin(const GeomVector<2, ownMemory>& lhs,
                                                                        const GeomVector<2, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<3> elementWiseMin(const GeomVector<3, ownMemory>& lhs,
                                                                        const GeomVector<3, otherMemory>& rhs);

template<bool ownMemory, bool otherMemory> GeomVector<1> elementWiseMax(const GeomVector<1, ownMemory>& lhs,
                                                                        const GeomVector<1, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<2> elementWiseMax(const GeomVector<2, ownMemory>& lhs,
                                                                        const GeomVector<2, otherMemory>& rhs);
template<bool ownMemory, bool otherMemory> GeomVector<3> elementWiseMax(const GeomVector<3, ownMemory>& lhs,
                                                                        const GeomVector<3, otherMemory>& rhs);

template<int nDim, bool ownMemory> std::istream& operator>>(std::istream& is, GeomVector<nDim, ownMemory>& vec);
template<int nDim, bool ownMemory> std::ostream& operator<<(std::ostream& os, const GeomVector<nDim, ownMemory>& vec);

}

#ifndef __GCCXML__
#include "GeomVectorInline.hh"
#endif

#else

// Forward declare the GeomVector class.
namespace Spheral {
  template<int nDim, bool ownMemory = true> class GeomVector;
}

#endif

