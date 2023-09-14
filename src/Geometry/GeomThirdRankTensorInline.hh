#include <algorithm>
#include <limits.h>
#include <string>

#include "GeomThirdRankTensor.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

#include "pybind11/pybind11.h"

//namespace PYBIND11_NAMESPACE {
//    template<int nDim> struct polymorphic_type_hook<Spheral::GeomThirdRankTensor<nDim>> {
//        static const void *get(const Spheral::GeomThirdRankTensor<nDim> *src, const std::type_info*& type) {
//            if (src) {
//                type = &typeid(Spheral::GeomThirdRankTensor<nDim>);
//                return static_cast<const Spheral::GeomThirdRankTensor<nDim>*>(src);
//            }
//            return src;
//        }
//    };
//} // namespace PYBIND11_NAMESPACE


namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
GeomThirdRankTensor():
  RankNTensor<nDim, 3, GeomThirdRankTensor>() {
}

//------------------------------------------------------------------------------
// Construct with the given value filling the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
GeomThirdRankTensor(const double val):
  RankNTensor<nDim, 3, GeomThirdRankTensor>(val) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
GeomThirdRankTensor(const GeomThirdRankTensor& rhs):
  RankNTensor<nDim, 3, GeomThirdRankTensor>(rhs) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
~GeomThirdRankTensor() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator=(const GeomThirdRankTensor& rhs) {
  RankNTensor<nDim, 3, GeomThirdRankTensor>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Assignment (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator=(const double rhs) {
  RankNTensor<nDim, 3, GeomThirdRankTensor>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomThirdRankTensor<nDim>::
operator()(const GeomThirdRankTensor::size_type i,
           const GeomThirdRankTensor::size_type j,
           const GeomThirdRankTensor::size_type k) const {
  REQUIRE(i < nDim and j < nDim and k < nDim);
  return mElements[(i*nDim + j)*nDim + k];
}

template<int nDim>
inline
double&
GeomThirdRankTensor<nDim>::
operator()(const GeomThirdRankTensor::size_type i,
           const GeomThirdRankTensor::size_type j,
           const GeomThirdRankTensor::size_type k) {
  REQUIRE(i < nDim and j < nDim and k < nDim);
  return mElements[(i*nDim + j)*nDim + k];
}

//------------------------------------------------------------------------------
// Have the global functions use the generic RankNTensor methods.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
operator*(const double lhs, const GeomThirdRankTensor<nDim>& rhs) {
  return lhs * dynamic_cast<const RankNTensor<nDim, 3, GeomThirdRankTensor<nDim> >&>(rhs);
}

template<int nDim>
inline
::std::istream&
operator>>(std::istream& is, GeomThirdRankTensor<nDim>& rhs) {
  return operator>>(is, dynamic_cast<RankNTensor<nDim, 3, GeomThirdRankTensor<nDim> >&>(rhs));
}

template<int nDim>
inline
::std::ostream&
operator<<(std::ostream& os, const GeomThirdRankTensor<nDim>& rhs) {
  return operator<<(os, dynamic_cast<const RankNTensor<nDim, 3, GeomThirdRankTensor<nDim> >&>(rhs));
}

}
