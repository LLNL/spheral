//---------------------------------Spheral++----------------------------------//
// GeomTensor -- Geometric Tensor Class.
//
// Created by JMO, Thu Apr 22 20:48:23 PDT 1999
// Modified by:
//   Aug 1, 99: JMO, removing MTL inheritance due to compilation problems,
//              and adding symmetric tensor type.
//   2004-08-19: JMO, converting to use boost::ublas library.
//   2004-08-23:  JMO, ublas is still too slow, so going to primitive C 
//                internal data types in accordance with suggestions from
//                Brian White
//   2004-11-10: JMO, experimenting with replacing the double array with individual
//               double elements, again for speed.
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeomTensor_default_hh__
#define __Spheral_GeomTensor_default_hh__

#include "Geometry/GeomVector_fwd.hh"
#include "Geometry/GeomTensor_fwd.hh"
#include "Geometry/GeomSymmetricTensor_fwd.hh"
#include "Geometry/GeomTensorBase.hh"

#include <iostream>
#include "Eigen/Dense"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace Spheral {

template<int nDim>
class GeomTensor: public GeomTensorBase<nDim> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;
  typedef Eigen::Matrix<double, nDim, nDim> EigenType;

  // Useful static memeber data.
  static const size_type nDimensions;
  static constexpr size_type numElements = nDim * nDim;
  static const GeomTensor zero;
  static const GeomTensor one;

  // Constructors.
  GeomTensor();
  explicit GeomTensor(const double a11);
  GeomTensor(const double a11, const double a12,
             const double a21, const double a22);
  GeomTensor(const double a11, const double a12, const double a13,
             const double a21, const double a22, const double a23,
             const double a31, const double a32, const double a33);
  GeomTensor(const GeomTensor& ten);
  explicit GeomTensor(const GeomSymmetricTensor<nDim>& ten);
  GeomTensor(const EigenType& ten);
  template<typename Derived> GeomTensor(const Eigen::MatrixBase<Derived>& ten);

  // Destructor.
  ~GeomTensor();

  // Assignment.
  GeomTensor& operator=(const GeomTensor& rhs);
  GeomTensor& operator=(const GeomSymmetricTensor<nDim>& rhs);
  template<typename Derived> GeomTensor& operator=(const Eigen::MatrixBase<Derived>& rhs);

  // Access the elements by indicies.
  double operator()(const size_type row, const size_type column) const;
  double& operator()(const size_type row, const size_type column);

  // More C++ style indexing.
  double operator[](size_type index) const;
  double& operator[](size_type index);

  // Access the individual elements by (x,y,z) notation.
  double xx() const;
  double xy() const;
  double xz() const;
  double yx() const;
  double yy() const;
  double yz() const;
  double zx() const;
  double zy() const;
  double zz() const;
  void xx(double val);
  void xy(double val);
  void xz(double val);
  void yx(double val);
  void yy(double val);
  void yz(double val);
  void zx(double val);
  void zy(double val);
  void zz(double val);

  // Get/set rows and columns.
  GeomVector<nDim> getRow(size_type index) const;
  GeomVector<nDim> getColumn(size_type index) const;
  void setRow(size_type index, const GeomVector<nDim>& vec);
  void setColumn(size_type index, const GeomVector<nDim>& vec);

  // Iterator access to the raw data.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // The zero and identity matrices.
  void Zero();
  void Identity();

  GeomTensor operator-() const;

  GeomTensor operator+(const GeomTensor& rhs) const;
  GeomTensor operator-(const GeomTensor& rhs) const;
  GeomTensor operator*(const GeomTensor& rhs) const;

  GeomTensor operator+(const GeomSymmetricTensor<nDim>& rhs) const;
  GeomTensor operator-(const GeomSymmetricTensor<nDim>& rhs) const;
  GeomTensor operator*(const GeomSymmetricTensor<nDim>& rhs) const;

  GeomVector<nDim> operator*(const GeomVector<nDim>& rhs) const;
  // GeomTensor operator+(const double val) const;
  // GeomTensor operator-(const double val) const;
  GeomTensor operator*(const double val) const;
  GeomTensor operator/(const double val) const;

  GeomTensor& operator+=(const GeomTensor& rhs);
  GeomTensor& operator-=(const GeomTensor& rhs);
  GeomTensor& operator*=(const GeomTensor& rhs);

  GeomTensor& operator+=(const GeomSymmetricTensor<nDim>& rhs);
  GeomTensor& operator-=(const GeomSymmetricTensor<nDim>& rhs);
  GeomTensor& operator*=(const GeomSymmetricTensor<nDim>& rhs);

  template<typename Derived> GeomTensor& operator+=(const Eigen::MatrixBase<Derived>& rhs);
  template<typename Derived> GeomTensor& operator-=(const Eigen::MatrixBase<Derived>& rhs);

  GeomTensor& operator*=(const double rhs);
  GeomTensor& operator/=(const double rhs);

  bool operator==(const GeomTensor& rhs) const;
  bool operator!=(const GeomTensor& rhs) const;
  bool operator<(const GeomTensor& rhs) const;
  bool operator>(const GeomTensor& rhs) const;
  bool operator<=(const GeomTensor& rhs) const;
  bool operator>=(const GeomTensor& rhs) const;

  bool operator==(const GeomSymmetricTensor<nDim>& rhs) const;
  bool operator!=(const GeomSymmetricTensor<nDim>& rhs) const;
  bool operator<(const GeomSymmetricTensor<nDim>& rhs) const;
  bool operator>(const GeomSymmetricTensor<nDim>& rhs) const;
  bool operator<=(const GeomSymmetricTensor<nDim>& rhs) const;
  bool operator>=(const GeomSymmetricTensor<nDim>& rhs) const;

  // bool operator==(const double rhs) const;
  // bool operator!=(const double rhs) const;
  // bool operator<(const double rhs) const;
  // bool operator>(const double rhs) const;
  // bool operator<=(const double rhs) const;
  // bool operator>=(const double rhs) const;

  GeomSymmetricTensor<nDim> Symmetric() const;
  GeomTensor SkewSymmetric() const;
  GeomTensor Transpose() const;
  GeomTensor Inverse() const;
  GeomVector<nDim> diagonalElements() const;
  double Trace() const;
  double Determinant() const;
  GeomVector<nDim> dot(const GeomVector<nDim>& vec) const;
  GeomTensor dot(const GeomTensor& rhs) const;
  GeomTensor dot(const GeomSymmetricTensor<nDim>& rhs) const;
  double doubledot(const GeomTensor& rhs) const;
  double doubledot(const GeomSymmetricTensor<nDim>& rhs) const;
  double selfDoubledot() const;

  // Return the square of this tensor (using matrix multiplication).
  GeomTensor square() const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  GeomTensor squareElements() const;

  // A simple method for returning the eigenvalues of a tensor.
  GeomVector<nDim> eigenValues() const;
  
  // Apply a rotational transform to this tensor ( R^-1 * (*this) * R ).
  void rotationalTransform(const GeomTensor& R);

  // Return the max absolute element.
  double maxAbsElement() const;

  //  Convert to an Eigen Vector
  EigenType eigen() const;

private:
  //--------------------------- Private Interface ---------------------------//
  size_type elementIndex(const size_type row, const size_type column) const;
};

// Declare specializations.
#ifndef WIN32
template<> const unsigned      GeomTensor<1>::nDimensions;
template<> const GeomTensor<1> GeomTensor<1>::zero;
template<> const GeomTensor<1> GeomTensor<1>::one;

template<> const unsigned      GeomTensor<2>::nDimensions;
template<> const GeomTensor<2> GeomTensor<2>::zero;
template<> const GeomTensor<2> GeomTensor<2>::one;

template<> const unsigned      GeomTensor<3>::nDimensions;
template<> const GeomTensor<3> GeomTensor<3>::zero;
template<> const GeomTensor<3> GeomTensor<3>::one;
#endif

template<> GeomTensor<2>::GeomTensor(const double, const double,
                                     const double, const double);
template<> GeomTensor<3>::GeomTensor(const double, const double, const double,
                                     const double, const double, const double,
                                     const double, const double, const double);

template<> GeomTensor<1>& GeomTensor<1>::operator=(const GeomTensor<1>& rhs);
template<> GeomTensor<2>& GeomTensor<2>::operator=(const GeomTensor<2>& rhs);
template<> GeomTensor<3>& GeomTensor<3>::operator=(const GeomTensor<3>& rhs);

template<> GeomTensor<1>& GeomTensor<1>::operator=(const GeomSymmetricTensor<1>& rhs);
template<> GeomTensor<2>& GeomTensor<2>::operator=(const GeomSymmetricTensor<2>& rhs);
template<> GeomTensor<3>& GeomTensor<3>::operator=(const GeomSymmetricTensor<3>& rhs);

template<> double GeomTensor<1>::xy() const;
template<> double GeomTensor<1>::xz() const;
template<> double GeomTensor<1>::yx() const;
template<> double GeomTensor<1>::yy() const;
template<> double GeomTensor<1>::yz() const;
template<> double GeomTensor<1>::zx() const;
template<> double GeomTensor<1>::zy() const;
template<> double GeomTensor<1>::zz() const;

template<> double GeomTensor<2>::xz() const;
template<> double GeomTensor<2>::yz() const;
template<> double GeomTensor<2>::zx() const;
template<> double GeomTensor<2>::zy() const;
template<> double GeomTensor<2>::zz() const;

template<> void GeomTensor<1>::xy(const double);
template<> void GeomTensor<1>::xz(const double);
template<> void GeomTensor<1>::yx(const double);
template<> void GeomTensor<1>::yy(const double);
template<> void GeomTensor<1>::yz(const double);
template<> void GeomTensor<1>::zx(const double);
template<> void GeomTensor<1>::zy(const double);
template<> void GeomTensor<1>::zz(const double);

template<> void GeomTensor<2>::xz(const double);
template<> void GeomTensor<2>::yz(const double);
template<> void GeomTensor<2>::zx(const double);
template<> void GeomTensor<2>::zy(const double);
template<> void GeomTensor<2>::zz(const double);

template<> GeomVector<1> GeomTensor<1>::getRow(const GeomTensor<1>::size_type) const;
template<> GeomVector<2> GeomTensor<2>::getRow(const GeomTensor<2>::size_type) const;
template<> GeomVector<3> GeomTensor<3>::getRow(const GeomTensor<3>::size_type) const;

template<> GeomVector<1> GeomTensor<1>::getColumn(const GeomTensor<1>::size_type) const;
template<> GeomVector<2> GeomTensor<2>::getColumn(const GeomTensor<2>::size_type) const;
template<> GeomVector<3> GeomTensor<3>::getColumn(const GeomTensor<3>::size_type) const;

template<> void GeomTensor<1>::setRow(const GeomTensor<1>::size_type, const GeomVector<1>&);
template<> void GeomTensor<2>::setRow(const GeomTensor<2>::size_type, const GeomVector<2>&);
template<> void GeomTensor<3>::setRow(const GeomTensor<3>::size_type, const GeomVector<3>&);

template<> void GeomTensor<1>::setColumn(const GeomTensor<1>::size_type, const GeomVector<1>&);
template<> void GeomTensor<2>::setColumn(const GeomTensor<2>::size_type, const GeomVector<2>&);
template<> void GeomTensor<3>::setColumn(const GeomTensor<3>::size_type, const GeomVector<3>&);

template<> GeomTensor<1> GeomTensor<1>::operator-() const;
template<> GeomTensor<2> GeomTensor<2>::operator-() const;
template<> GeomTensor<3> GeomTensor<3>::operator-() const;

// template<> GeomTensor<1> GeomTensor<1>::operator+(const double) const;
// template<> GeomTensor<2> GeomTensor<2>::operator+(const double) const;
// template<> GeomTensor<3> GeomTensor<3>::operator+(const double) const;

// template<> GeomTensor<1> GeomTensor<1>::operator-(const double) const;
// template<> GeomTensor<2> GeomTensor<2>::operator-(const double) const;
// template<> GeomTensor<3> GeomTensor<3>::operator-(const double) const;

template<> GeomTensor<1> GeomTensor<1>::operator*(const double) const;
template<> GeomTensor<2> GeomTensor<2>::operator*(const double) const;
template<> GeomTensor<3> GeomTensor<3>::operator*(const double) const;

template<> GeomTensor<1> GeomTensor<1>::operator/(const double) const;
template<> GeomTensor<2> GeomTensor<2>::operator/(const double) const;
template<> GeomTensor<3> GeomTensor<3>::operator/(const double) const;

template<> GeomTensor<1>& GeomTensor<1>::operator+=(const GeomTensor<1>&);
template<> GeomTensor<2>& GeomTensor<2>::operator+=(const GeomTensor<2>&);
template<> GeomTensor<3>& GeomTensor<3>::operator+=(const GeomTensor<3>&);

template<> GeomTensor<1>& GeomTensor<1>::operator+=(const GeomSymmetricTensor<1>&);
template<> GeomTensor<2>& GeomTensor<2>::operator+=(const GeomSymmetricTensor<2>&);
template<> GeomTensor<3>& GeomTensor<3>::operator+=(const GeomSymmetricTensor<3>&);

template<> GeomTensor<1>& GeomTensor<1>::operator-=(const GeomTensor<1>&);
template<> GeomTensor<2>& GeomTensor<2>::operator-=(const GeomTensor<2>&);
template<> GeomTensor<3>& GeomTensor<3>::operator-=(const GeomTensor<3>&);

template<> GeomTensor<1>& GeomTensor<1>::operator-=(const GeomSymmetricTensor<1>&);
template<> GeomTensor<2>& GeomTensor<2>::operator-=(const GeomSymmetricTensor<2>&);
template<> GeomTensor<3>& GeomTensor<3>::operator-=(const GeomSymmetricTensor<3>&);

template<> GeomTensor<1>& GeomTensor<1>::operator*=(const GeomTensor<1>&);
template<> GeomTensor<2>& GeomTensor<2>::operator*=(const GeomTensor<2>&);
template<> GeomTensor<3>& GeomTensor<3>::operator*=(const GeomTensor<3>&);
template<> GeomTensor<1>& GeomTensor<1>::operator*=(const GeomSymmetricTensor<1>&);
template<> GeomTensor<2>& GeomTensor<2>::operator*=(const GeomSymmetricTensor<2>&);
template<> GeomTensor<3>& GeomTensor<3>::operator*=(const GeomSymmetricTensor<3>&);

// template<> GeomTensor<1>& GeomTensor<1>::operator+=(const double);
// template<> GeomTensor<2>& GeomTensor<2>::operator+=(const double);
// template<> GeomTensor<3>& GeomTensor<3>::operator+=(const double);

// template<> GeomTensor<1>& GeomTensor<1>::operator-=(const double);
// template<> GeomTensor<2>& GeomTensor<2>::operator-=(const double);
// template<> GeomTensor<3>& GeomTensor<3>::operator-=(const double);

template<> GeomTensor<1>& GeomTensor<1>::operator*=(const double);
template<> GeomTensor<2>& GeomTensor<2>::operator*=(const double);
template<> GeomTensor<3>& GeomTensor<3>::operator*=(const double);

template<> GeomTensor<1>& GeomTensor<1>::operator/=(const double);
template<> GeomTensor<2>& GeomTensor<2>::operator/=(const double);
template<> GeomTensor<3>& GeomTensor<3>::operator/=(const double);

template<> bool GeomTensor<1>::operator==(const GeomTensor<1>&) const;
template<> bool GeomTensor<2>::operator==(const GeomTensor<2>&) const;
template<> bool GeomTensor<3>::operator==(const GeomTensor<3>&) const;

template<> bool GeomTensor<1>::operator==(const GeomSymmetricTensor<1>&) const;
template<> bool GeomTensor<2>::operator==(const GeomSymmetricTensor<2>&) const;
template<> bool GeomTensor<3>::operator==(const GeomSymmetricTensor<3>&) const;

template<> GeomSymmetricTensor<1> GeomTensor<1>::Symmetric() const;
template<> GeomSymmetricTensor<2> GeomTensor<2>::Symmetric() const;
template<> GeomSymmetricTensor<3> GeomTensor<3>::Symmetric() const;

template<> GeomTensor<1> GeomTensor<1>::SkewSymmetric() const;
template<> GeomTensor<2> GeomTensor<2>::SkewSymmetric() const;
template<> GeomTensor<3> GeomTensor<3>::SkewSymmetric() const;

template<> GeomTensor<1> GeomTensor<1>::Transpose() const;
template<> GeomTensor<2> GeomTensor<2>::Transpose() const;
template<> GeomTensor<3> GeomTensor<3>::Transpose() const;

template<> GeomTensor<1> GeomTensor<1>::Inverse() const;
template<> GeomTensor<2> GeomTensor<2>::Inverse() const;
template<> GeomTensor<3> GeomTensor<3>::Inverse() const;

template<> GeomVector<1> GeomTensor<1>::diagonalElements() const;
template<> GeomVector<2> GeomTensor<2>::diagonalElements() const;
template<> GeomVector<3> GeomTensor<3>::diagonalElements() const;

template<> double GeomTensor<1>::Trace() const;
template<> double GeomTensor<2>::Trace() const;
template<> double GeomTensor<3>::Trace() const;

template<> double GeomTensor<1>::Determinant() const;
template<> double GeomTensor<2>::Determinant() const;
template<> double GeomTensor<3>::Determinant() const;

template<> GeomVector<1> GeomTensor<1>::dot(const GeomVector<1>&) const;
template<> GeomVector<2> GeomTensor<2>::dot(const GeomVector<2>&) const;
template<> GeomVector<3> GeomTensor<3>::dot(const GeomVector<3>&) const;

template<> GeomTensor<1> GeomTensor<1>::dot(const GeomTensor<1>&) const;
template<> GeomTensor<2> GeomTensor<2>::dot(const GeomTensor<2>&) const;
template<> GeomTensor<3> GeomTensor<3>::dot(const GeomTensor<3>&) const;
template<> GeomTensor<1> GeomTensor<1>::dot(const GeomSymmetricTensor<1>&) const;
template<> GeomTensor<2> GeomTensor<2>::dot(const GeomSymmetricTensor<2>&) const;
template<> GeomTensor<3> GeomTensor<3>::dot(const GeomSymmetricTensor<3>&) const;

template<> double GeomTensor<1>::doubledot(const GeomTensor<1>&) const;
template<> double GeomTensor<2>::doubledot(const GeomTensor<2>&) const;
template<> double GeomTensor<3>::doubledot(const GeomTensor<3>&) const;

template<> double GeomTensor<1>::doubledot(const GeomSymmetricTensor<1>&) const;
template<> double GeomTensor<2>::doubledot(const GeomSymmetricTensor<2>&) const;
template<> double GeomTensor<3>::doubledot(const GeomSymmetricTensor<3>&) const;

template<> GeomTensor<1> GeomTensor<1>::square() const;
template<> GeomTensor<2> GeomTensor<2>::square() const;
template<> GeomTensor<3> GeomTensor<3>::square() const;

template<> GeomTensor<1> GeomTensor<1>::squareElements() const;
template<> GeomTensor<2> GeomTensor<2>::squareElements() const;
template<> GeomTensor<3> GeomTensor<3>::squareElements() const;

template<> void GeomTensor<1>::rotationalTransform(const GeomTensor<1>&);
template<> void GeomTensor<2>::rotationalTransform(const GeomTensor<2>&);
template<> void GeomTensor<3>::rotationalTransform(const GeomTensor<3>&);

template<> double GeomTensor<1>::maxAbsElement() const;
template<> double GeomTensor<2>::maxAbsElement() const;
template<> double GeomTensor<3>::maxAbsElement() const;

template<> const GeomTensor<1> GeomTensor<1>::zero;
template<> const GeomTensor<2> GeomTensor<2>::zero;
template<> const GeomTensor<3> GeomTensor<3>::zero;

// Forward declare the global functions.
template<int nDim> GeomTensor<nDim> operator*(double lhs, const GeomTensor<nDim>& rhs);
template<int nDim> ::std::istream& operator>>(std::istream& is, GeomTensor<nDim>& ten);
template<int nDim> std::ostream& operator<<(std::ostream& os, const GeomTensor<nDim>& ten);

#if defined(_OPENMP) && _OPENMP >= 201107
#pragma omp declare reduction(tensadd : GeomTensor<1> : omp_out += omp_in ) initializer( omp_priv = GeomTensor<1>(0.0) )
#pragma omp declare reduction(tensdif : GeomTensor<1> : omp_out -= omp_in ) initializer( omp_priv = GeomTensor<1>(0.0) )
#pragma omp declare reduction(tensadd : GeomTensor<2> : omp_out += omp_in ) initializer( omp_priv = GeomTensor<2>(0.0,0.0,0.0,0.0)) 
#pragma omp declare reduction(tensdif : GeomTensor<2> : omp_out -= omp_in ) initializer( omp_priv = GeomTensor<2>(0.0,0.0,0.0,0.0))
#pragma omp declare reduction(tensadd : GeomTensor<3> : omp_out += omp_in ) initializer( omp_priv = GeomTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) )
#pragma omp declare reduction(tensdif : GeomTensor<3> : omp_out -= omp_in ) initializer( omp_priv = GeomTensor<3>(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0) ) 
#endif

}

#include "GeomTensorInline_default.hh"

#endif
