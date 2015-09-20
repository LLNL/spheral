//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#include "TableKernel.hh"

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/bisectSearch.hh"

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

inline
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

inline
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

inline
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
// Construct from a kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename KernelType>
TableKernel<Dimension>::TableKernel(const KernelType& kernel,
                                    const int numPoints,
                                    const double hmult):
  Kernel<Dimension, TableKernel<Dimension> >(),
  mKernelValues(0),
  mGradValues(0),
  mGrad2Values(0),
  mAkernel(0),
  mBkernel(0),
  mAgrad(0),
  mBgrad(0),
  mAgrad2(0),
  mBgrad2(0),
  mNumPoints(0),
  mStepSize(0.0),
  mNperhValues(0),
  mWsumValues(0),
  mMinNperh(0.5),
  mMaxNperh(10.0) {

  // Pre-conditions.
  VERIFY(numPoints > 0);
  VERIFY(hmult > 0.0);

  // Set the volume normalization and kernel extent.
  this->setVolumeNormalization(kernel.volumeNormalization() / Dimension::pownu(hmult));
  this->setKernelExtent(hmult * kernel.kernelExtent());
  this->setInflectionPoint(hmult * kernel.inflectionPoint());

  // Set the number of points and table step size.
  CHECK(numPoints > 1);
  mNumPoints = numPoints;
  mStepSize = this->kernelExtent()/(numPoints - 1);
  CHECK(stepSize() > 0.0);

  // Resize the tables for our data.
  mKernelValues = std::vector<double>(numPoints);
  mGradValues = std::vector<double>(numPoints);
  mGrad2Values = std::vector<double>(numPoints);
  mAkernel = std::vector<double>(numPoints);
  mBkernel = std::vector<double>(numPoints);
  mAgrad = std::vector<double>(numPoints);
  mBgrad = std::vector<double>(numPoints);
  mAgrad2 = std::vector<double>(numPoints);
  mBgrad2 = std::vector<double>(numPoints);

  // Fill in the kernel and gradient values.  Note that we will go ahead and fold
  // the normalization constants in here, so we don't have to multiply by them
  // in the value lookups.
  const double deta = mStepSize/hmult;
  for (int i = 0; i < numPoints; ++i) {
    CHECK(i*mStepSize >= 0.0);
    mKernelValues[i] = kernel(i*deta, 1.0);
    mGradValues[i] = kernel.grad(i*deta, 1.0);
    mGrad2Values[i] = kernel.grad2(i*deta, 1.0);
  }

  // Set the delta kernel values for internal use.
  this->setParabolicCoeffs();

  // // Adjust the table kernel values to reproduce the volume integral as closely as possible.
  // const double vol1 = simpsonsVolumeIntegral<Dimension, TableKernel<Dimension> >(*this, 0.0, kernel.kernelExtent(), 10*numPoints);
  // const double f = 1.0/vol1; // /(this->volumeNormalization());
  // std::cerr << "Scaling: " << f << " " << vol1 << " " << this->volumeNormalization() << std::endl;
  // for (int i = 0; i < numPoints; ++i) {
  //   mKernelValues[i] *= f;
  // }

  // // Reset the delta kernel values for internal use.
  // this->setDeltaKernelValues();

  // Set the table of n per h values.
  this->setNperhValues();

  // That should be it, so we should have left the kernel in a valid state.
  ENSURE(valid());
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
    mAkernel = rhs.mAkernel;
    mBkernel = rhs.mBkernel;
    mAgrad = rhs.mAgrad;
    mBgrad = rhs.mBgrad;
    mAgrad2 = rhs.mAgrad2;
    mBgrad2 = rhs.mBgrad2;
    mNumPoints = rhs.mNumPoints;
    mStepSize = rhs.mStepSize;
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
        (ub == mNumPoints and Wsum >= mWsumValues[mNumPoints - 1]) ||
        (Wsum >= mWsumValues[lb] and Wsum <= mWsumValues[ub]));

  // Now interpolate for the corresponding nodes per h (within bounds);
  double result;
  if (lb == -1) {
    result = mNperhValues[0];
  } else if (ub == mNumPoints) {
    result = mNperhValues[mNumPoints - 1];
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
        (ub == mNumPoints and nPerh >= mNperhValues[mNumPoints - 1]) ||
        (nPerh >= mNperhValues[lb] and nPerh <= mNperhValues[ub]));

  // Now interpolate for the corresponding Wsum.
  double result;
  if (lb == -1) {
    result = mWsumValues[0];
  } else if (ub == mNumPoints) {
    result = mWsumValues[mNumPoints - 1];
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
// Initialize the parabolic fit coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TableKernel<Dimension>::
setParabolicCoeffs() {
  
  // Size stuff up.
  REQUIRE(mNumPoints > 0);
  mAkernel = std::vector<double>(mNumPoints);
  mBkernel = std::vector<double>(mNumPoints);
  mAgrad = std::vector<double>(mNumPoints);
  mBgrad = std::vector<double>(mNumPoints);
  mAgrad2 = std::vector<double>(mNumPoints);
  mBgrad2 = std::vector<double>(mNumPoints);

  // Find the coefficient fits.
  for (int i0 = 0; i0 < mNumPoints - 2; ++i0) {
    const int i1 = i0 + 1;
    const int i2 = i0 + 2;
    CHECK(i2 < mNumPoints);

    mAkernel[i1] = ((mKernelValues[i2] - mKernelValues[i1])/mStepSize - (mKernelValues[i1] - mKernelValues[i0])/mStepSize)/(2.0*mStepSize);
    mBkernel[i1] = ((mKernelValues[i2] - mKernelValues[i1]) + (mKernelValues[i1] - mKernelValues[i0]))/(2.0*mStepSize);

    mAgrad[i1] = ((mGradValues[i2] - mGradValues[i1])/mStepSize - (mGradValues[i1] - mGradValues[i0])/mStepSize)/(2.0*mStepSize);
    mBgrad[i1] = ((mGradValues[i2] - mGradValues[i1]) + (mGradValues[i1] - mGradValues[i0]))/(2.0*mStepSize);

    mAgrad2[i1] = ((mGrad2Values[i2] - mGrad2Values[i1])/mStepSize - (mGrad2Values[i1] - mGrad2Values[i0])/mStepSize)/(2.0*mStepSize);
    mBgrad2[i1] = ((mGrad2Values[i2] - mGrad2Values[i1]) + (mGrad2Values[i1] - mGrad2Values[i0]))/(2.0*mStepSize);
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
  REQUIRE(mNumPoints > 1);
  REQUIRE(this->kernelExtent() > 0.0);

  // Size the Nperh array.
  mWsumValues.resize(mNumPoints);
  mNperhValues.resize(mNumPoints);

  // For the allowed range of n per h, sum up the kernel values.
  const double dnperh = (mMaxNperh - mMinNperh)/(mNumPoints - 1);
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
  for (int i = 0; i != mNumPoints - 1; ++i) {
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
          mKernelValues.size() == mNumPoints and
          mGradValues.size() == mNumPoints and
          mGrad2Values.size() == mNumPoints and
          mAkernel.size() == mNumPoints and
          mBkernel.size() == mNumPoints and
          mAgrad.size() == mNumPoints and
          mBgrad.size() == mNumPoints and
          mAgrad2.size() == mNumPoints and
          mBgrad2.size() == mNumPoints and
          mStepSize > 0.0);
}

}
}

