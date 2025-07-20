#include <algorithm>
#include <limits.h>
#include <string>

#include "GeomFourthRankTensor.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomFourthRankTensor<nDim>::
GeomFourthRankTensor():
  RankNTensor<nDim, 4, GeomFourthRankTensor>() {
}

//------------------------------------------------------------------------------
// Construct with the given value filling the tensor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomFourthRankTensor<nDim>::
GeomFourthRankTensor(const double val):
  RankNTensor<nDim, 4, GeomFourthRankTensor>(val) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomFourthRankTensor<nDim>::
GeomFourthRankTensor(const GeomFourthRankTensor& rhs):
  RankNTensor<nDim, 4, GeomFourthRankTensor>(rhs) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomFourthRankTensor<nDim>::
~GeomFourthRankTensor() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomFourthRankTensor<nDim>&
GeomFourthRankTensor<nDim>::
operator=(const GeomFourthRankTensor& rhs) {
  RankNTensor<nDim, 4, GeomFourthRankTensor>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Assignment (double).
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
GeomFourthRankTensor<nDim>&
GeomFourthRankTensor<nDim>::
operator=(const double rhs) {
  RankNTensor<nDim, 4, GeomFourthRankTensor>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
SPHERAL_HOST_DEVICE
inline
double
GeomFourthRankTensor<nDim>::
operator()(const GeomFourthRankTensor::size_type i,
           const GeomFourthRankTensor::size_type j,
           const GeomFourthRankTensor::size_type k,
           const GeomFourthRankTensor::size_type m) const {
  REQUIRE(i < nDim and j < nDim and k < nDim and m < nDim);
  return mElements[((i*nDim + j)*nDim + k)*nDim + m];
}

template<int nDim>
SPHERAL_HOST_DEVICE
inline
double&
GeomFourthRankTensor<nDim>::
operator()(const GeomFourthRankTensor::size_type i,
           const GeomFourthRankTensor::size_type j,
           const GeomFourthRankTensor::size_type k,
           const GeomFourthRankTensor::size_type m) {
  REQUIRE(i < nDim and j < nDim and k < nDim and m < nDim);
  return mElements[((i*nDim + j)*nDim + k)*nDim + m];
}

//------------------------------------------------------------------------------
// Have the global functions use the generic RankNTensor methods.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomFourthRankTensor<nDim>
operator*(const double lhs, const GeomFourthRankTensor<nDim>& rhs) {
  return lhs * dynamic_cast<const RankNTensor<nDim, 4, GeomFourthRankTensor<nDim> >&>(rhs);
}

template<int nDim>
inline
::std::istream&
operator>>(std::istream& is, GeomFourthRankTensor<nDim>& rhs) {
  return operator>>(is, dynamic_cast<const RankNTensor<nDim, 4, GeomFourthRankTensor<nDim> >&>(rhs));
}

template<int nDim>
inline
::std::ostream&
operator<<(std::ostream& os, const GeomFourthRankTensor<nDim>& rhs) {
  return operator<<(os, dynamic_cast<const RankNTensor<nDim, 4, GeomFourthRankTensor<nDim> >&>(rhs));
}

}
