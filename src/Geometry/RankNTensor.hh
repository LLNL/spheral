//---------------------------------Spheral++----------------------------------//
// RankNTensor -- Arbitrary rank tensor.  Assumes each nth index spans 
//                Dimension::nDim.
//
// Created by JMO, Sun Oct 11 10:38:26 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_RankNTensor_hh__
#define __Spheral_RankNTensor_hh__

#include <iostream>
#include "Utilities/FastMath.hh"

namespace Spheral {

template<int nDim, int rank, typename Descendant>
class RankNTensor {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef const double* const_iterator;
  typedef double* iterator;
  typedef unsigned size_type;

  // Useful static member data.
  static const size_type nrank;
  static const size_type nDimensions;
  static constexpr size_type numElements = FastMath::calcPower(nDim, rank);

  // Constructors.
  RankNTensor();
  explicit RankNTensor(const double val);
  RankNTensor(const RankNTensor& rhs);

  // Destructor.
  virtual ~RankNTensor();

  // Assignment.
  RankNTensor& operator=(const RankNTensor& rhs);
  RankNTensor& operator=(const double rhs);

  // More C++ style indexing.
  double operator[](size_type index) const;
  double& operator[](size_type index);

  // Iterator access to the raw data.
  iterator begin();
  iterator end();

  const_iterator begin() const;
  const_iterator end() const;

  // Zero out the tensor.
  void Zero();

  // Assorted operations.
  Descendant operator-() const;

  Descendant& operator+=(const RankNTensor& rhs);
  Descendant& operator-=(const RankNTensor& rhs);

  Descendant operator+(const RankNTensor& rhs) const;
  Descendant operator-(const RankNTensor& rhs) const;

  Descendant& operator*=(const double rhs);
  Descendant& operator/=(const double rhs);

  Descendant operator*(const double rhs) const;
  Descendant operator/(const double rhs) const;

  bool operator==(const RankNTensor& rhs) const;
  bool operator!=(const RankNTensor& rhs) const;
  bool operator<(const RankNTensor& rhs) const;
  bool operator>(const RankNTensor& rhs) const;
  bool operator<=(const RankNTensor& rhs) const;
  bool operator>=(const RankNTensor& rhs) const;

  bool operator==(const double rhs) const;
  bool operator!=(const double rhs) const;

  double doubledot(const RankNTensor<nDim, rank, Descendant>& rhs) const;

  // Return a tensor where each element is the square of the corresponding 
  // element of this tensor.
  Descendant squareElements() const;

  // Return the max absolute element.
  double maxAbsElement() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  double mElements[numElements] = {};
};

// Forward declare the global functions.
template<int nDim, int rank, typename Descendant> Descendant operator*(const double lhs, const RankNTensor<nDim, rank, Descendant>& rhs);

template<int nDim, int rank, typename Descendant> ::std::istream& operator>>(std::istream& is, RankNTensor<nDim, rank, Descendant>& ten);
template<int nDim, int rank, typename Descendant> ::std::ostream& operator<<(std::ostream& os, const RankNTensor<nDim, rank, Descendant>& ten);

// Initialize our static variables.
template<int nDim, int rank, typename Descendant> const typename RankNTensor<nDim, rank, Descendant>::size_type RankNTensor<nDim, rank, Descendant>::nrank = rank;
template<int nDim, int rank, typename Descendant> const typename RankNTensor<nDim, rank, Descendant>::size_type RankNTensor<nDim, rank, Descendant>::nDimensions = nDim;

}

#ifndef __GCCXML__
#include "RankNTensorInline.hh"
#endif

#endif
