//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#include "TableKernel.hh"

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/bisectSearch.hh"

#include "BSplineKernel.hh"
#include "W4SplineKernel.hh"
#include "GaussianKernel.hh"
#include "SuperGaussianKernel.hh"
#include "PiGaussianKernel.hh"
#include "SincKernel.hh"
#include "NSincPolynomialKernel.hh"
#include "NBSplineKernel.hh"
#include "HatKernel.hh"
#include "QuarticSplineKernel.hh"
#include "QuinticSplineKernel.hh"

namespace Spheral {
namespace KernelSpace {

using namespace std;

//------------------------------------------------------------------------------
// Sum the Kernel values for the given stepsize.
//------------------------------------------------------------------------------
// double
// sumKernelValues(const TableKernel<Dim<1> >& W,
//                 const double deta) {
//   REQUIRE(deta > 0);
//   double result = 0.0;
//   double etax = deta;
//   while (etax < W.kernelExtent()) {
//     result += 2.0*W(etax, 1.0);
//     etax += deta;
//   }
//   return result;
// }

// double
// sumKernelValues(const TableKernel<Dim<2> >& W,
//                 const double deta) {
//   REQUIRE(deta > 0);
//   double result = 0.0;
//   double etax = deta;
//   while (etax < W.kernelExtent()) {
//     result += 2.0*W(etax, 1.0);
//     etax += deta;
//   }
//   return 2.0*M_PI*result;
// }

// double
// sumKernelValues(const TableKernel<Dim<3> >& W,
//                 const double deta) {
//   REQUIRE(deta > 0);
//   double result = 0.0;
//   double etax = deta;
//   while (etax < W.kernelExtent()) {
//     result += 2.0*W(etax, 1.0);
//     etax += deta;
//   }
//   return 4.0*M_PI*result;
// }

double
sumKernelValues(const TableKernel<Dim<1> >& W,
                const double deta) {
  REQUIRE(deta > 0);
  double result = 0.0;
  double etax = deta;
  while (etax < W.kernelExtent()) {
    result += 2.0*std::abs(W.gradValue(etax, 1.0));
    etax += deta;
  }
  return result;
}

double
sumKernelValues(const TableKernel<Dim<2> >& W,
                const double deta) {
  REQUIRE(deta > 0);
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  double result = 0.0;
  double etay = 0.0;
  while (etay < W.kernelExtent()) {
    double etax = 0.0;
    while (etax < W.kernelExtent()) {
      const Vector eta(etax, etay);
      double dresult = std::abs(W.gradValue(eta.magnitude(), 1.0));
      if (distinctlyGreaterThan(etax, 0.0)) dresult *= 2.0;
      if (distinctlyGreaterThan(etay, 0.0)) dresult *= 2.0;
      if (fuzzyEqual(eta.magnitude(), 0.0)) dresult *= 0.0;
      result += dresult;
      etax += deta;
    }
    etay += deta;
  }
  return sqrt(result);
}

double
sumKernelValues(const TableKernel<Dim<3> >& W,
                const double deta) {
  REQUIRE(deta > 0);
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  double result = 0.0;
  double etaz = 0.0;
  while (etaz < W.kernelExtent()) {
    double etay = 0.0;
    while (etay < W.kernelExtent()) {
      double etax = 0.0;
      while (etax < W.kernelExtent()) {
        const Vector eta(etax, etay, etaz);
        CHECK(eta >= 0.0);
        double dresult = std::abs(W.gradValue(eta.magnitude(), 1.0));
        if (distinctlyGreaterThan(etax, 0.0)) dresult *= 2.0;
        if (distinctlyGreaterThan(etay, 0.0)) dresult *= 2.0;
        if (distinctlyGreaterThan(etaz, 0.0)) dresult *= 2.0;
        if (fuzzyEqual(eta.magnitude(), 0.0)) dresult *= 0.0;
        result += dresult;
        etax += deta;
      }
      etay += deta;
    }
    etaz += deta;
  }
  return FastMath::CubeRootHalley2(result);
}

// Special hacked version to allow for running 1-D stacks of nodes in 3-D.
// This is way ugly and tricky -- DON'T EMULATE THIS KIND OF EXAMPLE!
template<typename KernelType>
inline
double
sumKernelValuesAs1D(const KernelType& W,
                    const double deta) {
  REQUIRE(deta > 0);
  double result = 0.0;
  double etax = deta;
  while (etax < W.kernelExtent()) {
    result += 2.0*std::abs(W.gradValue(etax, 1.0));
    etax += deta;
  }
  return result;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>::
~TableKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
TableKernel<Dimension>&
TableKernel<Dimension>::
operator=(const TableKernel<Dimension>& rhs) {
  if (this != &rhs) {
    Kernel<Dimension, TableKernel<Dimension> >::operator=(rhs);
    mKernelValues = rhs.mKernelValues;
    mGradValues = rhs.mGradValues;
    mGrad2Values = rhs.mGrad2Values;
    mDeltaKernelValues = rhs.mDeltaKernelValues;
    mDeltaGradValues = rhs.mDeltaGradValues;
    mDeltaGrad2Values = rhs.mDeltaGrad2Values;
    mNumPoints = rhs.mNumPoints;
    mNumPoints1 = rhs.mNumPoints1;
    mStepSizeInv = rhs.mStepSizeInv;
    mNperhValues = rhs.mNperhValues;
    mWsumValues =  rhs.mWsumValues;
    mMinNperh = rhs.mMinNperh;
    mMaxNperh = rhs.mMaxNperh;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Determine the number of nodes per smoothing scale implied by the given
// sum of kernel values.
//------------------------------------------------------------------------------
template<typename Dimension>
double
TableKernel<Dimension>::
equivalentNodesPerSmoothingScale(const double Wsum) const {

  // Find the lower bound in the tabulated Wsum's bracketing the input
  // value.
  const int lb = bisectSearch(mWsumValues, Wsum);
  CHECK((lb >= -1) and (lb <= int(mWsumValues.size()) - 1));
  const int ub = lb + 1;
  CHECK((lb == -1 and Wsum <= mWsumValues[0]) ||
        (ub == mNumPoints and Wsum >= mWsumValues[mNumPoints1]) ||
        (Wsum >= mWsumValues[lb] and Wsum <= mWsumValues[ub]));

  // Now interpolate for the corresponding nodes per h (within bounds);
  double result;
  if (lb == -1) {
    result = mNperhValues[0];
  } else if (ub == mNumPoints) {
    result = mNperhValues[mNumPoints1];
  } else {
    result = std::min(mNperhValues[ub],
                      std::max(mNperhValues[lb],
                               mNperhValues[lb] +
                               (Wsum - mWsumValues[lb])/
                               (mWsumValues[ub] - mWsumValues[lb])*
                               (mNperhValues[ub] - mNperhValues[lb])));
    ENSURE(result >= mNperhValues[lb] and result <= mNperhValues[ub]);
  }
  return result;
}

//------------------------------------------------------------------------------
// Determine the effective Wsum we would expect for the given n per h.
//------------------------------------------------------------------------------
template<typename Dimension>
double
TableKernel<Dimension>::
equivalentWsum(const double nPerh) const {

  // Find the lower bound in the tabulated n per h's bracketing the input
  // value.
  const int lb = bisectSearch(mNperhValues, nPerh);
  CHECK((lb >= -1) and (lb <= int(mNperhValues.size()) - 1));
  const int ub = lb + 1;
  CHECK((lb == -1 and nPerh <= mNperhValues[0]) ||
        (ub == mNumPoints and nPerh >= mNperhValues[mNumPoints1]) ||
        (nPerh >= mNperhValues[lb] and nPerh <= mNperhValues[ub]));

  // Now interpolate for the corresponding Wsum.
  double result;
  if (lb == -1) {
    result = mWsumValues[0];
  } else if (ub == mNumPoints) {
    result = mWsumValues[mNumPoints1];
  } else {
    result = std::min(mWsumValues[ub], 
                      std::max(mWsumValues[lb],
                               mWsumValues[lb] +
                               (nPerh - mNperhValues[lb])/
                               (mNperhValues[ub] - mNperhValues[lb])*
                               (mWsumValues[ub] - mWsumValues[lb])));
    ENSURE(result >= mWsumValues[lb] and result <= mWsumValues[ub]);
  }
  return result;
}

//------------------------------------------------------------------------------
// Initialize the delta kernel values.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TableKernel<Dimension>::
setDeltaKernelValues() {

  // Size the delta fields.
  REQUIRE(mNumPoints1 >= 0 and
          mNumPoints1 == mNumPoints - 1);
  mDeltaKernelValues.resize(mNumPoints1);
  mDeltaGradValues.resize(mNumPoints1);
  mDeltaGrad2Values.resize(mNumPoints1);

  // Now fill in the deltas.
  CHECK(mKernelValues.size() == mNumPoints);
  CHECK(mGradValues.size() == mNumPoints);
  CHECK(mGrad2Values.size() == mNumPoints);
  for (int i = 0; i < mNumPoints1; ++i) {
    int j = i + 1;
    mDeltaKernelValues[i] = mKernelValues[j] - mKernelValues[i];
    mDeltaGradValues[i] = mGradValues[j] - mGradValues[i];
    mDeltaGrad2Values[i] = mGrad2Values[j] - mGrad2Values[i];
  }
}

//------------------------------------------------------------------------------
// Initialize the Nperh values.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TableKernel<Dimension>::
setNperhValues(const bool scaleTo1D) {
  REQUIRE(mMinNperh > 0.0);
  REQUIRE(mMaxNperh > mMinNperh);
  REQUIRE(mNumPoints1 > 0);
  REQUIRE(this->kernelExtent() > 0.0);

  // Size the Nperh array.
  mWsumValues.resize(mNumPoints);
  mNperhValues.resize(mNumPoints);

  // For the allowed range of n per h, sum up the kernel values.
  const double dnperh = (mMaxNperh - mMinNperh)/mNumPoints1;
  for (int i = 0; i != mNumPoints; ++i) {
    const double nperh = mMinNperh + i*dnperh;
    CHECK(nperh >= mMinNperh and nperh <= mMaxNperh);
    const double deta = 1.0/nperh;
    mNperhValues[i] = nperh;
    if (scaleTo1D) {
      mWsumValues[i] = sumKernelValuesAs1D(*this, deta);
    } else {
      mWsumValues[i] = sumKernelValues(*this, deta);
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  ENSURE(mWsumValues.size() == mNumPoints);
  ENSURE(mNperhValues.size() == mNumPoints);
  for (int i = 0; i != mNumPoints1; ++i) {
    ENSURE(mWsumValues[i] <= mWsumValues[i + 1]);
    ENSURE(mNperhValues[i] <= mNperhValues[i + 1]);
  }
  END_CONTRACT_SCOPE;

}

//------------------------------------------------------------------------------
// Determine if the kernel is in a valid, ready to use state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TableKernel<Dimension>::
valid() const {
  return (Kernel<Dimension, TableKernel<Dimension> >::valid() and
          mNumPoints > 0 and
          mNumPoints1 == mNumPoints - 1 and
          mKernelValues.size() == mNumPoints and
          mGradValues.size() == mNumPoints and
          mGrad2Values.size() == mNumPoints and
          mDeltaKernelValues.size() == mNumPoints - 1 and
          mDeltaGradValues.size() == mNumPoints - 1 and
          mDeltaGrad2Values.size() == mNumPoints - 1 and
          mStepSizeInv > 0.0);
}

}
}

