#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"
#include "VolumeIntegrationFunctions.hh"

namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
TableKernel<Dimension>::kernelValue(const double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  return Hdet*parabolicInterp(etaMagnitude, mKernelValues, mAkernel, mBkernel);
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
  return Hdet*parabolicInterp(etaMagnitude, mGradValues, mAgrad, mBgrad);
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
  return Hdet*parabolicInterp(etaMagnitude, mGrad2Values, mAgrad2, mBgrad2);
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
  return std::make_pair(Hdet*parabolicInterp(etaMagnitude, mKernelValues, mAkernel, mBkernel),
                        Hdet*parabolicInterp(etaMagnitude, mGradValues, mAgrad, mBgrad));
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
    kernelValues[i] = Hdets[i]*parabolicInterp(etaMagnitudes[i], mKernelValues, mAkernel, mBkernel);
    gradValues[i] = Hdets[i]*parabolicInterp(etaMagnitudes[i], mGradValues, mAgrad, mBgrad);
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
  return parabolicInterp(etaMagnitude, mf1Values, mAf1, mBf1);
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
  return parabolicInterp(etaMagnitude, mf2Values, mAf2, mBf2);
}

//------------------------------------------------------------------------------
// Return the f1 and f2 RZ values for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::pair<double, double>
TableKernel<Dimension>::f1Andf2(const double etaMagnitude) const {
  VERIFY2(false, "TableKernel::f1Andf2 lookup only valid for 2D kernels.");
}

template<>
inline
std::pair<double, double>
TableKernel<Dim<2> >::f1Andf2(const double etaMagnitude) const {
  REQUIRE(etaMagnitude >= 0.0);
  return std::make_pair(parabolicInterp(etaMagnitude, mf1Values, mAf1, mBf2),
                        parabolicInterp(etaMagnitude, mf2Values, mAf2, mBf2));
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
                                        const std::vector<double>& table,
                                        const std::vector<double>& a,
                                        const std::vector<double>& b) const {
  REQUIRE(table.size() == mNumPoints);
  REQUIRE(a.size() == mNumPoints);
  REQUIRE(b.size() == mNumPoints);
  if (etaMagnitude < this->mKernelExtent) {
    const int i0 = min(mNumPoints - 3, lowerBound(etaMagnitude));
    const int i1 = i0 + 1;
    const double deta = etaMagnitude - i1*mStepSize;
    return a[i1]*deta*deta + b[i1]*deta + table[i1];
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the assorted tabular lookup data.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
kernelValues() const {
  return mKernelValues;
}

template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
gradValues() const {
  return mGradValues;
}

template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
grad2Values() const {
  return mGrad2Values;
}

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

template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
f1Values() const {
  return mf1Values;
}

template<typename Dimension>
inline
const std::vector<double>&
TableKernel<Dimension>::
f2Values() const {
  return mf2Values;
}

}
}

