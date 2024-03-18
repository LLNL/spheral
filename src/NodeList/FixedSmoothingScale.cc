//---------------------------------Spheral++----------------------------------//
// FixedSmoothingScale
//
// Implements the static fixed smoothing scale option.
//
// Created by JMO, Wed Sep 14 13:50:49 PDT 2005
//----------------------------------------------------------------------------//
#include "FixedSmoothingScale.hh"
#include "Field/FieldList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FixedSmoothingScale<Dimension>::
FixedSmoothingScale():
  SmoothingScaleBase<Dimension>() {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FixedSmoothingScale<Dimension>::
FixedSmoothingScale(const FixedSmoothingScale<Dimension>& rhs):
  SmoothingScaleBase<Dimension>(rhs) {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
FixedSmoothingScale<Dimension>&
FixedSmoothingScale<Dimension>::
operator=(const FixedSmoothingScale& rhs) {
  SmoothingScaleBase<Dimension>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
FixedSmoothingScale<Dimension>::
~FixedSmoothingScale<Dimension>() {
}

//------------------------------------------------------------------------------
// Time derivative of the smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
FixedSmoothingScale<Dimension>::
smoothingScaleDerivative(const SymTensor& /*H*/,
                         const Vector& /*pos*/,
                         const Tensor& /*DvDx*/,
                         const Scalar /*hmin*/,
                         const Scalar /*hmax*/,
                         const Scalar /*hminratio*/,
                         const Scalar /*nPerh*/) const {
  return SymTensor::zero;
}

//------------------------------------------------------------------------------
// Directly evaluate the smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
FixedSmoothingScale<Dimension>::
newSmoothingScale(const SymTensor& H,
                  const Vector& /*pos*/,
                  const Scalar /*zerothMoment*/,
                  const Vector& /*firstMoment*/,
                  const SymTensor& /*secondMomentEta*/,
                  const SymTensor& /*secondMomentLab*/,
                  const TableKernel<Dimension>& /*W*/,
                  const Scalar /*hmin*/,
                  const Scalar /*hmax*/,
                  const Scalar /*hminratio*/,
                  const Scalar /*nPerh*/,
                  const ConnectivityMap<Dimension>& /*connectivityMap*/,
                  const unsigned /*nodeListi*/,
                  const unsigned /*i*/) const {
  return H;
}

//------------------------------------------------------------------------------
// Directly evaluate the smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
FixedSmoothingScale<Dimension>::
idealSmoothingScale(const SymTensor& H,
                    const Vector& /*pos*/,
                    const Scalar /*zerothMoment*/,
                    const Vector& /*firstMoment*/,
                    const SymTensor& /*secondMomentEta*/,
                    const SymTensor& /*secondMomentLab*/,
                    const TableKernel<Dimension>& /*W*/,
                    const Scalar /*hmin*/,
                    const Scalar /*hmax*/,
                    const Scalar /*hminratio*/,
                    const Scalar /*nPerh*/,
                    const ConnectivityMap<Dimension>& /*connectivityMap*/,
                    const unsigned /*nodeListi*/,
                    const unsigned /*i*/) const {
  return H;
}

//------------------------------------------------------------------------------
// Use the volumes of tessellation to set the new Hs.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::SymTensor
FixedSmoothingScale<Dimension>::
idealSmoothingScale(const SymTensor& H,
                    const Mesh<Dimension>& /*mesh*/,
                    const typename Mesh<Dimension>::Zone& /*zone*/,
                    const Scalar /*hmin*/,
                    const Scalar /*hmax*/,
                    const Scalar /*hminratio*/,
                    const Scalar /*nPerh*/) const {
  return H;
}

}
