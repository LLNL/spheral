#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "VolumeIntegrationFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::kernelValue(const double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaMagnitude < this->mKernelExtent) {
    return Hdet*parabolicInterp(etaMagnitude, mAkernel, mBkernel, mCkernel);
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
TableKernel<Dimension>::gradValue(const double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaMagnitude < this->mKernelExtent) {
    return Hdet*parabolicInterp(etaMagnitude, mAgrad, mBgrad, mCgrad);
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
TableKernel<Dimension>::grad2Value(const double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaMagnitude < this->mKernelExtent) {
    return Hdet*parabolicInterp(etaMagnitude, mAgrad2, mBgrad2, mCgrad2);
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
TableKernel<Dimension>::kernelAndGradValue(const double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaMagnitude < this->mKernelExtent) {
    const int i0 = std::min(mNumPoints - 3, lowerBound(etaMagnitude));
    const int i1 = i0 + 1;
    CHECK(i1 >= 1 and i1 <= mNumPoints - 2);
    const double x = etaMagnitude/mStepSize - i0;
    CHECK(x >= 0.0);
    return std::make_pair(Hdet*(mAkernel[i1] + mBkernel[i1]*x + mCkernel[i1]*x*x),
                          Hdet*(mAgrad[i1] + mBgrad[i1]*x + mCgrad[i1]*x*x));
  } else {
    return std::make_pair(0.0, 0.0);
  }
  // return std::make_pair(Hdet*parabolicInterp(etaMagnitude, mAkernel, mBkernel, mCkernel),
  //                       Hdet*parabolicInterp(etaMagnitude, mAgrad, mBgrad, mCgrad));
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
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(Hdets.size() == n);
    for (size_t i = 0; i != n; ++i) {
      REQUIRE(etaMagnitudes[i] >= 0.0);
      REQUIRE(Hdets[i] >= 0.0);
    }
  }
  END_CONTRACT_SCOPE

  // Prepare the results.
  kernelValues = std::vector<double>(n);
  gradValues = std::vector<double>(n);

  // Fill those suckers in.
  for (size_t i = 0; i != n; ++i) {
    kernelValues[i] = Hdets[i]*parabolicInterp(etaMagnitudes[i], mAkernel, mBkernel, mCkernel);
    gradValues[i] = Hdets[i]*parabolicInterp(etaMagnitudes[i], mAgrad, mBgrad, mCgrad);
  }
}

//------------------------------------------------------------------------------
// Return the f1 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::f1(const double etaMagnitude) const {
  VERIFY2(false, "TableKernel::f1 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::f1(const double etaMagnitude) const {
  REQUIRE(etaMagnitude >= 0.0);
  if (etaMagnitude < this->mKernelExtent - mStepSize) {
    return parabolicInterp(etaMagnitude, mAf1, mBf1, mCf1);
  } else {
    return 1.0;
  }
}

//------------------------------------------------------------------------------
// Return the f2 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::f2(const double etaMagnitude) const {
  VERIFY2(false, "TableKernel::f2 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::f2(const double etaMagnitude) const {
  REQUIRE(etaMagnitude >= 0.0);
  return parabolicInterp(etaMagnitude, mAf2, mBf2, mCf2);
}

//------------------------------------------------------------------------------
// Return the grad_r f1 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::gradf1(const double etaMagnitude) const {
  VERIFY2(false, "TableKernel::gradf1 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::gradf1(const double etaMagnitude) const {
  REQUIRE(etaMagnitude >= 0.0);
  return parabolicInterp(etaMagnitude, mAgradf1, mBgradf1, mCgradf1);
}

//------------------------------------------------------------------------------
// Return the grad_r f2 RZ correctiong for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::gradf2(const double etaMagnitude) const {
  VERIFY2(false, "TableKernel::gradf2 lookup only valid for 2D kernels.");
}

template<>
inline
double
TableKernel<Dim<2> >::gradf2(const double etaMagnitude) const {
  REQUIRE(etaMagnitude >= 0.0);
  return parabolicInterp(etaMagnitude, mAgradf2, mBgradf2, mCgradf2);
}

//------------------------------------------------------------------------------
// Return the (f1, f2, gradf1, gradf2) RZ values for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
TableKernel<Dimension>::f1Andf2(const double etaMagnitude,
                                double& f1,
                                double& f2,
                                double& gradf1,
                                double& gradf2) const {
  VERIFY2(false, "TableKernel::f1Andf2 lookup only valid for 2D kernels.");
}

template<>
inline
void
TableKernel<Dim<2> >::f1Andf2(const double etaMagnitude,
                              double& f1,
                              double& f2,
                              double& gradf1,
                              double& gradf2) const {
  REQUIRE(etaMagnitude >= 0.0);
  if (etaMagnitude < this->mKernelExtent - mStepSize) {
    const int i0 = std::min(mNumPoints - 3, lowerBound(etaMagnitude));
    const int i1 = i0 + 1;
    CHECK(i1 >= 1 and i1 <= mNumPoints - 2);
    const double x = etaMagnitude/mStepSize - i0;
    CHECK(x >= 0.0);
    f1 = mAf1[i1] + mBf1[i1]*x + mCf1[i1]*x*x;
    f2 = mAf2[i1] + mBf2[i1]*x + mCf2[i1]*x*x;
    gradf1 = mAgradf1[i1] + mBgradf1[i1]*x + mCgradf1[i1]*x*x;
    gradf2 = mAgradf2[i1] + mBgradf2[i1]*x + mCgradf2[i1]*x*x;
  } else {
    f1 = 1.0;
    f2 = 1.0;
    gradf1 = 0.0;
    gradf2 = 0.0;
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
  return mStepSize;
}

template<typename Dimension>
inline
double
TableKernel<Dimension>::stepSizeInv() const {
  REQUIRE(mStepSize > 0.0);
  return 1.0/mStepSize;
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
  const int result = std::min(mNumPoints - 1, int(etaMagnitude/mStepSize));
  ENSURE(result >= 0 && result < mNumPoints);
  return result;
}

//------------------------------------------------------------------------------
// Helper to do the parabolic interpolation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::parabolicInterp(const double etaMagnitude, 
                                        const std::vector<double>& a,
                                        const std::vector<double>& b,
                                        const std::vector<double>& c) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(a.size() == mNumPoints);
  REQUIRE(b.size() == mNumPoints);
  REQUIRE(c.size() == mNumPoints);
  const int i0 = std::min(mNumPoints - 3, lowerBound(etaMagnitude));
  const int i1 = i0 + 1;
  CHECK(i1 >= 1 and i1 <= mNumPoints - 2);
  const double x = etaMagnitude/mStepSize - i0;
  CHECK(x >= 0.0);
  return a[i1] + b[i1]*x + c[i1]*x*x;
}

//------------------------------------------------------------------------------
// Return the assorted tabular lookup data.
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

