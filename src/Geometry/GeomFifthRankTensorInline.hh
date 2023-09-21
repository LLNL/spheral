#include <algorithm>
#include <limits.h>
#include <string>

#include "GeomFifthRankTensor.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFifthRankTensor<nDim>::
GeomFifthRankTensor():
  RankNTensor<nDim, 5, GeomFifthRankTensor>() {
}

//------------------------------------------------------------------------------
// Construct with the given value filling the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFifthRankTensor<nDim>::
GeomFifthRankTensor(const double val):
  RankNTensor<nDim, 5, GeomFifthRankTensor>(val) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFifthRankTensor<nDim>::
GeomFifthRankTensor(const GeomFifthRankTensor& rhs):
  RankNTensor<nDim, 5, GeomFifthRankTensor>(rhs) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFifthRankTensor<nDim>::
~GeomFifthRankTensor() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFifthRankTensor<nDim>&
GeomFifthRankTensor<nDim>::
operator=(const GeomFifthRankTensor& rhs) {
  RankNTensor<nDim, 5, GeomFifthRankTensor>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Assignment (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFifthRankTensor<nDim>&
GeomFifthRankTensor<nDim>::
operator=(const double rhs) {
  RankNTensor<nDim, 5, GeomFifthRankTensor>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomFifthRankTensor<nDim>::
operator()(const GeomFifthRankTensor::size_type i,
           const GeomFifthRankTensor::size_type j,
           const GeomFifthRankTensor::size_type k,
           const GeomFifthRankTensor::size_type m,
           const GeomFifthRankTensor::size_type n) const {
  REQUIRE(i < nDim and j < nDim and k < nDim and m < nDim and n < nDim);
  return mElements[(((i*nDim + j)*nDim + k)*nDim + m)*nDim + n];
}

template<int nDim>
inline
double&
GeomFifthRankTensor<nDim>::
operator()(const GeomFifthRankTensor::size_type i,
           const GeomFifthRankTensor::size_type j,
           const GeomFifthRankTensor::size_type k,
           const GeomFifthRankTensor::size_type m,
           const GeomFifthRankTensor::size_type n) {
  REQUIRE(i < nDim and j < nDim and k < nDim and m < nDim and n < nDim);
  return mElements[(((i*nDim + j)*nDim + k)*nDim + m)*nDim + n];
}

//------------------------------------------------------------------------------
// Have the global functions use the generic RankNTensor methods.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFifthRankTensor<nDim>
operator*(const double lhs, const GeomFifthRankTensor<nDim>& rhs) {
  return lhs * dynamic_cast<const RankNTensor<nDim, 5, GeomFifthRankTensor<nDim> >&>(rhs);
}

template<int nDim>
inline
::std::istream&
operator>>(std::istream& is, GeomFifthRankTensor<nDim>& rhs) {
  return operator>>(is, dynamic_cast<const RankNTensor<nDim, 5, GeomFifthRankTensor<nDim> >&>(rhs));
}

template<int nDim>
inline
::std::ostream&
operator<<(std::ostream& os, const GeomFifthRankTensor<nDim>& rhs) {
  return operator<<(os, dynamic_cast<const RankNTensor<nDim, 5, GeomFifthRankTensor<nDim> >&>(rhs));
}

}
