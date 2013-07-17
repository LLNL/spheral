//---------------------------------Spheral++----------------------------------//
// GeomThirdRankTensor -- Geometric Tensor (rank 3) Class.
//
// This is a very simple, limited functionality rank 3 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim.
//
// Created by JMO, Thu Jul 17 17:05:58 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomThirdRankTensor_hh__
#define __Spheral_GeomThirdRankTensor_hh__

#include <iostream>

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"
#include "Geometry/GeomThirdRankTensor_fwd.hh"

namespace Spheral {

template<int nDim>
class GeomThirdRankTensor {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;

  // Useful static member data.
  static const size_type nDimensions;
  static const size_type numElements;
  static const GeomThirdRankTensor zero;

  // Constructors.
  GeomThirdRankTensor();
  explicit GeomThirdRankTensor(const double val);
  GeomThirdRankTensor(const GeomThirdRankTensor& ten);

  // Destructor.
  ~GeomThirdRankTensor();

  // Assignment.
  GeomThirdRankTensor& operator=(const GeomThirdRankTensor& rhs);
  GeomThirdRankTensor& operator=(const double rhs);

  // Access the elements by indicies.
  double operator()(const size_type i, const size_type j, const size_type k) const;
  double& operator()(const size_type i, const size_type j, const size_type k);

  // Iterator access to the raw data.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // Zero out the tensor.
  void Zero();

  // Assorted operations.
  GeomThirdRankTensor operator-() const;

  GeomThirdRankTensor& operator+=(const GeomThirdRankTensor& rhs);
  GeomThirdRankTensor& operator-=(const GeomThirdRankTensor& rhs);

  GeomThirdRankTensor operator+(const GeomThirdRankTensor& rhs) const;
  GeomThirdRankTensor operator-(const GeomThirdRankTensor& rhs) const;

  GeomThirdRankTensor& operator*=(const double rhs);
  GeomThirdRankTensor& operator/=(const double rhs);

  GeomThirdRankTensor operator*(const double rhs) const;
  GeomThirdRankTensor operator/(const double rhs) const;

  bool operator==(const GeomThirdRankTensor& rhs) const;
  bool operator!=(const GeomThirdRankTensor& rhs) const;
  bool operator<(const GeomThirdRankTensor& rhs) const;
  bool operator>(const GeomThirdRankTensor& rhs) const;
  bool operator<=(const GeomThirdRankTensor& rhs) const;
  bool operator>=(const GeomThirdRankTensor& rhs) const;

  bool operator==(const double rhs) const;
  bool operator!=(const double rhs) const;

  double doubledot(const GeomThirdRankTensor& rhs) const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  GeomThirdRankTensor squareElements() const;

  // Return the max absolute element.
  double maxAbsElement() const;

private:
  //--------------------------- Private Interface ---------------------------//
  size_type elementIndex(const size_type i,
                         const size_type j,
                         const size_type k) const;

  double* mElements;
  static const size_type nDim2;
};

// Forward declare the global functions.
template<int nDim> GeomThirdRankTensor<nDim> operator*(const double lhs, const GeomThirdRankTensor<nDim>& rhs);

template<int nDim> ::std::istream& operator>>(std::istream& is, GeomThirdRankTensor<nDim>& ten);
template<int nDim> ::std::ostream& operator<<(std::ostream& os, const GeomThirdRankTensor<nDim>& ten);

}

#ifndef __GCCXML__
#include "GeomThirdRankTensorInline.hh"
#endif

#endif
