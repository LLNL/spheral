#include <math.h>
#include "Utilities/DBC.hh"
#include "Geometry/GeomVector.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return as a reference descendent class.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
Descendent&
Kernel<Dimension, Descendent>::asDescendent() const {
  return static_cast<Descendent&>(const_cast<Kernel<Dimension, Descendent>&>(*this));
//   return static_cast<Descendent&>(*this);
}

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
Kernel<Dimension, Descendent>::Kernel():
  mVolumeNormalization(0.0),
  mKernelExtent(0.0),
  mInflectionPoint(0.0) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
Kernel<Dimension, Descendent>::~Kernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
Kernel<Dimension, Descendent>&
Kernel<Dimension, Descendent>::operator=(const Kernel<Dimension, Descendent>& rhs) {
  if (this != &rhs) {
    mVolumeNormalization = rhs.volumeNormalization();
    mKernelExtent = rhs.kernelExtent();
    mInflectionPoint = rhs.inflectionPoint();
  }
  return *this;
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance or position.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::operator()(double etaMagnitude, 
                                          const HType& H) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return kernelValue(etaMagnitude, H.Determinant());
}

template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::operator()(const typename Dimension::Vector& eta,
                                          const HType& H) const {
  REQUIRE(valid());
  return kernelValue(eta.magnitude(), H.Determinant());
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::operator()(double etaMagnitude, 
                                          const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return kernelValue(etaMagnitude, Hdet);
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::operator()(const typename Dimension::Vector& eta,
                                          const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  return kernelValue(eta.magnitude(), Hdet);
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::grad(double etaMagnitude, const HType& H) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return gradValue(etaMagnitude, H.Determinant());
}

template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::grad(const typename Dimension::Vector& eta,
                                    const HType& H) const {
  REQUIRE(valid());
  return gradValue(eta.magnitude(), H.Determinant());
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::grad(double etaMagnitude,
                                    const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return gradValue(etaMagnitude, Hdet);
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::grad(const typename Dimension::Vector& eta,
                                    const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  return gradValue(eta.magnitude(), Hdet);
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::grad2(double etaMagnitude, const HType& H) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return grad2Value(etaMagnitude, H.Determinant());
}

template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::grad2(const typename Dimension::Vector& eta,
                                     const HType& H) const {
  REQUIRE(valid());
  return grad2Value(eta.magnitude(), H.Determinant());
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::grad2(double etaMagnitude,
                                     const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return grad2Value(etaMagnitude, Hdet);
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::grad2(const typename Dimension::Vector& eta,
                                     const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  return grad2Value(eta.magnitude(), Hdet);
}

//------------------------------------------------------------------------------
// Return the gradient with respect to h value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::gradh(double etaMagnitude, const HType& H) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return gradhValue(etaMagnitude, H.Determinant());
}

template<typename Dimension, typename Descendent>
template<typename HType>
inline
double
Kernel<Dimension, Descendent>::gradh(const typename Dimension::Vector& eta,
                                     const HType& H) const {
  REQUIRE(valid());
  return gradhValue(eta.magnitude(), H.Determinant());
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::gradh(double etaMagnitude,
                                     const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  REQUIRE(etaMagnitude >= 0.0);
  return gradhValue(etaMagnitude, Hdet);
}

template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::gradh(const typename Dimension::Vector& eta,
                                     const typename Dimension::Scalar& Hdet) const {
  REQUIRE(valid());
  return gradhValue(eta.magnitude(), Hdet);
}

//------------------------------------------------------------------------------
// All kernels must redefine this method to give the W value for a given
// distance and H determinant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::kernelValue(double etaMagnitude, double Hdet) const {
  REQUIRE(valid());
  return asDescendent().kernelValue(etaMagnitude, Hdet);
}

//------------------------------------------------------------------------------
// All kernels must redefine this method to give the gradient value for a given
// distance and H determinant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::gradValue(double etaMagnitude, double Hdet) const {
  REQUIRE(valid());
  return asDescendent().gradValue(etaMagnitude, Hdet);
}

//------------------------------------------------------------------------------
// All kernels must redefine this method to give the second derivative value for
// a given distance and H determinant.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::grad2Value(double etaMagnitude, double Hdet) const {
  REQUIRE(valid());
  return asDescendent().grad2Value(etaMagnitude, Hdet);
}

//------------------------------------------------------------------------------
// Compute the gradient with respect to h.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::gradhValue(double etaMagnitude, double Hdet) const {
  REQUIRE(valid());
  return -etaMagnitude * Dimension::rootnu(Hdet) * asDescendent().gradValue(etaMagnitude, Hdet);
}

//------------------------------------------------------------------------------
// Return the volume normalization.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::volumeNormalization() const {
  return mVolumeNormalization;
}

//------------------------------------------------------------------------------
// Set the volume normalization.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
void
Kernel<Dimension, Descendent>::setVolumeNormalization(double volumeNormalization) {
  REQUIRE(volumeNormalization > 0.0);
  mVolumeNormalization = volumeNormalization;
}

//------------------------------------------------------------------------------
// Return the kernel extent.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::kernelExtent() const {
  return mKernelExtent;
}

//------------------------------------------------------------------------------
// Set the kernel extent.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
void
Kernel<Dimension, Descendent>::setKernelExtent(double extent) {
  REQUIRE(extent > 0.0);
  mKernelExtent = extent;
}

//------------------------------------------------------------------------------
// Return the inflection point.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
double
Kernel<Dimension, Descendent>::inflectionPoint() const {
  return mInflectionPoint;
}

//------------------------------------------------------------------------------
// Set the inflection point.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
void
Kernel<Dimension, Descendent>::setInflectionPoint(double x) {
  mInflectionPoint = x;
}

//------------------------------------------------------------------------------
// Test if the kernel is valid.
//------------------------------------------------------------------------------
template<typename Dimension, typename Descendent>
inline
bool
Kernel<Dimension, Descendent>::valid() const {
  return (volumeNormalization() > 0.0 && kernelExtent() > 0.0);
}

}
