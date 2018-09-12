#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the order of the density fit in a cell.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
VoronoiHourglassControl<Dimension>::
order() const {
  return mOrder;
}

template<typename Dimension>
inline
void
VoronoiHourglassControl<Dimension>::
order(const unsigned x) {
  VERIFY(x >= 0 and x <= 1);
  mOrder = x;
}

//------------------------------------------------------------------------------
// Access the slope limiter for the density gradient in a cell.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
VoronoiHourglassControl<Dimension>::
limiter() const {
  return mLimiter;
}

template<typename Dimension>
inline
void
VoronoiHourglassControl<Dimension>::
limiter(const unsigned x) {
  VERIFY(x >= 0 and x <= 2);
  mLimiter = x;
}

//------------------------------------------------------------------------------
// Access the fraction of relaxation we allow.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
VoronoiHourglassControl<Dimension>::
fraction() const {
  return mFraction;
}

template<typename Dimension>
inline
void
VoronoiHourglassControl<Dimension>::
fraction(const double x) {
  VERIFY(x >= 0);
  mFraction = x;
}

//------------------------------------------------------------------------------
// Access the mask determining which nodes we filter.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
VoronoiHourglassControl<Dimension>::
mask() const {
  return mMask;
}

template<typename Dimension>
inline
void
VoronoiHourglassControl<Dimension>::
mask(const FieldList<Dimension, int>& x) {
  VERIFY(x.localMin() >= 0 and x.localMax() <= 1);
  mMask = x;
}

//------------------------------------------------------------------------------
// Access the kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
VoronoiHourglassControl<Dimension>::
kernel() const {
  return mW;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
VoronoiHourglassControl<Dimension>::
gradRho() const {
  return mGradRho;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
VoronoiHourglassControl<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
VoronoiHourglassControl<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
VoronoiHourglassControl<Dimension>::
C() const {
  return mC;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
VoronoiHourglassControl<Dimension>::
D() const {
  return mD;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
VoronoiHourglassControl<Dimension>::
gradA() const {
  return mGradA;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Tensor>&
VoronoiHourglassControl<Dimension>::
gradB() const {
  return mGradB;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
VoronoiHourglassControl<Dimension>::
weight() const {
  return mWeight;
}

}
