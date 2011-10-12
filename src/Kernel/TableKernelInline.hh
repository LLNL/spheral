#include "Geometry/Dimension.hh"
#include "DBC.hh"

namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
// Construct from a kernel.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename KernelType>
inline
TableKernel<Dimension>::TableKernel(const KernelType& kernel,
                                    const int numPoints,
                                    const double hmult):
  Kernel<Dimension, TableKernel<Dimension> >(),
  mKernelValues(0),
  mGradValues(0),
  mGrad2Values(0),
  mDeltaKernelValues(0),
  mDeltaGradValues(0),
  mDeltaGrad2Values(0),
  mNumPoints(0),
  mNumPoints1(-1),
  mStepSizeInv(0.0),
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
  mNumPoints1 = mNumPoints - 1;
  mStepSizeInv = 1.0/(this->kernelExtent()/(numPoints - 1));
  CHECK(stepSize() > 0.0);

  // Resize the tables for our data.
  mKernelValues.resize(numPoints);
  mGradValues.resize(numPoints);
  mGrad2Values.resize(numPoints);
  CHECK(mKernelValues.size() == numPoints);
  CHECK(mGradValues.size() == numPoints);
  CHECK(mGrad2Values.size() == numPoints);

  // Fill in the kernel and gradient values.  Note that we will go ahead and fold
  // the normalization constants in here, so we don't have to multiply by them
  // in the value lookups.
  const double deta = stepSize() / hmult;
  for (int i = 0; i < numPoints; ++i) {
    CHECK(i*stepSize() >= 0.0);
    mKernelValues[i] = kernel(i*deta, 1.0);
    mGradValues[i] = kernel.grad(i*deta, 1.0);
    mGrad2Values[i] = kernel.grad2(i*deta, 1.0);
  }

  // Set the delta kernel values for internal use.
  this->setDeltaKernelValues();

  // Set the table of n per h values.
  this->setNperhValues();

  // That should be it, so we should have left the kernel in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::kernelValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < this->mKernelExtent) {
    const int iLower = lowerBound(etaMagnitude);
    const int iUpper = std::min(iLower + 1, mNumPoints1);
    const double thpt = std::max(0.0, std::min(1.0, etaMagnitude*stepSizeInv() - iLower));
    CHECK(iUpper >= iLower);
    CHECK(thpt >= 0.0 and thpt <= 1.0);
    return Hdet*(mKernelValues[iLower] + mDeltaKernelValues[iLower]*thpt);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::gradValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < this->mKernelExtent) {
    const int iLower = lowerBound(etaMagnitude);
    const int iUpper = std::min(iLower + 1, mNumPoints1);
    const double thpt = std::max(0.0, std::min(1.0, etaMagnitude*stepSizeInv() - iLower));
    CHECK(iUpper >= iLower);
    CHECK(thpt >= 0.0 and thpt <= 1.0);
    return Hdet*(mGradValues[iLower] + mDeltaGradValues[iLower]*thpt);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::grad2Value(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < this->mKernelExtent) {
    const int iLower = lowerBound(etaMagnitude);
    const int iUpper = std::min(iLower + 1, mNumPoints1);
    const double thpt = std::max(0.0, std::min(1.0, etaMagnitude*stepSizeInv() - iLower));
    CHECK(iUpper >= iLower);
    CHECK(thpt >= 0.0 and thpt <= 1.0);
    return Hdet*(mGrad2Values[iLower] + mDeltaGrad2Values[iLower]*thpt);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::pair<double, double>
TableKernel<Dimension>::kernelAndGradValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < this->mKernelExtent) {
    const int iLower = lowerBound(etaMagnitude);
    const int iUpper = std::min(iLower + 1, mNumPoints1);
    const double thpt = std::max(0.0, std::min(1.0, etaMagnitude*stepSizeInv() - iLower));
    CHECK(iUpper >= iLower);
    CHECK(thpt >= 0.0 and thpt <= 1.0);
    return std::pair<double, double>
      (Hdet*(mKernelValues[iLower] + mDeltaKernelValues[iLower]*thpt),
       Hdet*(mGradValues[iLower] + mDeltaGradValues[iLower]*thpt));
  } else {
    return std::pair<double, double>(0.0, 0.0);
  }
}

//------------------------------------------------------------------------------
// Return the kernel and gradient values for a set of normalized distances.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TableKernel<Dimension>::kernelAndGradValues(const std::vector<double>& etaMagnitudes,
                                            const std::vector<double>& Hdets,
                                            std::vector<double>& kernelValues,
                                            std::vector<double>& gradValues) const {
  // Preconditions.
  const size_t n = etaMagnitudes.size();
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(Hdets.size() == n);
    for (size_t i = 0; i != n; ++i) {
      REQUIRE(etaMagnitudes[i] >= 0.0);
      REQUIRE(Hdets[i] >= 0.0);
    }
  }
  END_CONTRACT_SCOPE;

  // Prepare the results.
  kernelValues.resize(n);
  gradValues.resize(n);

  // Fill those suckers in.
  for (size_t i = 0; i != n; ++i) {
    const double& etaMagnitude = etaMagnitudes[i];
    const double& Hdet = Hdets[i];
    if (etaMagnitude < this->mKernelExtent) {
      const int iLower = lowerBound(etaMagnitude);
      const int iUpper = std::min(iLower + 1, mNumPoints1);
      const double thpt = std::max(0.0, std::min(1.0, etaMagnitude*stepSizeInv() - iLower));
      CHECK(iUpper >= iLower);
      CHECK(thpt >= 0.0 and thpt <= 1.0);
      kernelValues[i] = Hdet*(mKernelValues[iLower] + mDeltaKernelValues[iLower]*thpt);
      gradValues[i] = Hdet*(mGradValues[iLower] + mDeltaGradValues[iLower]*thpt);
    } else {
      kernelValues[i] = 0.0;
      gradValues[i] = 0.0;
    }
  }
}

//------------------------------------------------------------------------------
// Return the number of points in the tables.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
TableKernel<Dimension>::numPoints() const {
  return mNumPoints;
}

//------------------------------------------------------------------------------
// Return the table step size in eta.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::stepSize() const {
  REQUIRE(mStepSizeInv > 0.0);
  return 1.0/mStepSizeInv;
}

template<typename Dimension>
inline
double
TableKernel<Dimension>::stepSizeInv() const {
  return mStepSizeInv;
}

//------------------------------------------------------------------------------
// Return the lower bound index in the table for the given normalized radial
// position.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
TableKernel<Dimension>::lowerBound(double etaMagnitude) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(stepSizeInv() > 0.0);
  const int result = std::min(mNumPoints1, int(etaMagnitude*mStepSizeInv));
  ENSURE(result >= 0 && result < numPoints());
  return result;
}

//------------------------------------------------------------------------------
// Return the tabular lookup data for the nodes per smoothing scale.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
nperhValues() const {
  return mNperhValues;
}

template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
WsumValues() const {
  return mWsumValues;
}

}
}

