//---------------------------------Spheral++----------------------------------//
// TableKernel -- Build an interpolation kernel using interpolation between
// tabulated points.
//
// Created by JMO, Mon Jun 19 21:06:28 PDT 2000
//----------------------------------------------------------------------------//
#include "Eigen/Dense"

#include "TableKernel.hh"

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/safeInv.hh"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Sum the Kernel values for the given stepsize.
//------------------------------------------------------------------------------
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
// Compute the (f1,f2) integrals relation for the given zeta = r/h 
// (RZ corrections).
//------------------------------------------------------------------------------
template<typename KernelType>
class volfunc {
  const KernelType& W;
public:
  volfunc(const KernelType& W): W(W) {}
  double operator()(const double eta) const {
    return W.kernelValue(eta, 1.0);
  }
};

template<typename KernelType>
class f1func {
  const KernelType& W;
  double zeta;
public:
  f1func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    return abs(safeInvVar(zeta)*eta)*W.kernelValue(abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class f2func {
  const KernelType& W;
  double zeta;
public:
  f2func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    return safeInvVar(zeta*zeta)*eta*abs(eta)*W.kernelValue(abs(zeta - eta), 1.0);
  }
};


template<typename KernelType>
class gradf1func {
  const KernelType& W;
  double zeta;
public:
  gradf1func(const KernelType& W, const double zeta): W(W), zeta(zeta) {}
  double operator()(const double eta) const {
    const double Wu = W.kernelValue(abs(zeta - eta), 1.0);
    const double gWu = W.gradValue(abs(zeta - eta), 1.0);
    const double gf1inv = safeInvVar(zeta)*abs(eta)*gWu - safeInvVar(zeta*zeta)*abs(eta)*Wu;
    if (eta < 0.0) {
      return -gf1inv;
    } else {
      return gf1inv;
    }
  }
};


template<typename KernelType>
double
f1Integral(const KernelType& W,
           const double zeta,
           const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return safeInvVar(simpsonsIntegration<f1func<KernelType>, double, double>(f1func<KernelType>(W, zeta), 
                                                                            zeta - etaMax, 
                                                                            zeta + etaMax,
                                                                            numbins));
}

template<typename KernelType>
double
f2Integral(const KernelType& W,
           const double zeta,
           const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return safeInvVar(simpsonsIntegration<f2func<KernelType>, double, double>(f2func<KernelType>(W, zeta), 
                                                                            zeta - etaMax, 
                                                                            zeta + etaMax,
                                                                            numbins));
}

template<typename KernelType>
double
gradf1Integral(const KernelType& W,
               const double zeta,
               const unsigned numbins) {
  const double etaMax = W.kernelExtent();
  CHECK(zeta <= etaMax);
  return simpsonsIntegration<gradf1func<KernelType>, double, double>(gradf1func<KernelType>(W, zeta), 
                                                                     zeta - etaMax, 
                                                                     zeta + etaMax,
                                                                     numbins);
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
  mAkernel(),
  mBkernel(),
  mCkernel(),
  mAgrad(),
  mBgrad(),
  mCgrad(),
  mAgrad2(),
  mBgrad2(),
  mCgrad2(),
  mNumPoints(0),
  mStepSize(0.0),
  mNperhValues(),
  mWsumValues(),
  mMinNperh(0.25),
  mMaxNperh(10.0),
  mAf1(),
  mBf1(),
  mCf1(),
  mAf2(),
  mBf2(),
  mCf2() {

  // Pre-conditions.
  VERIFY(numPoints > 0);
  VERIFY(hmult > 0.0);

  // Set the volume normalization and kernel extent.
  this->setVolumeNormalization(1.0); // (kernel.volumeNormalization() / Dimension::pownu(hmult));  // We now build this into the tabular kernel values.
  this->setKernelExtent(hmult * kernel.kernelExtent());
  this->setInflectionPoint(hmult * kernel.inflectionPoint());

  // Set the number of points and table step size.
  CHECK(numPoints > 1);
  mNumPoints = numPoints;
  mStepSize = this->kernelExtent()/(numPoints - 1);
  CHECK(stepSize() > 0.0);

  // Set the parabolic fit coefficients for kernel and it's gradient.
  // Note that we will go ahead and fold the normalization constants in here, 
  // so we don't have to multiply by them in the value lookups.
  std::vector<double> kernelValues(numPoints);
  std::vector<double> gradValues(numPoints);
  std::vector<double> grad2Values(numPoints);
  const double correction = 1.0/Dimension::pownu(hmult);
  const double deta = mStepSize/hmult;
  for (int i = 0; i < numPoints; ++i) {
    CHECK(i*mStepSize >= 0.0);
    kernelValues[i] = correction*kernel(i*deta, 1.0);
    gradValues[i] = correction*kernel.grad(i*deta, 1.0);
    grad2Values[i] = correction*kernel.grad2(i*deta, 1.0);
  }
  setParabolicCoeffs(kernelValues, mAkernel, mBkernel, mCkernel);
  setParabolicCoeffs(gradValues, mAgrad, mBgrad, mCgrad);
  setParabolicCoeffs(grad2Values, mAgrad2, mBgrad2, mCgrad2);

  // If we're a 2D kernel we set the RZ correction information.
  if (Dimension::nDim == 2) {
    std::vector<double> f1Values(numPoints);
    std::vector<double> f2Values(numPoints);
    std::vector<double> gradf1Values(numPoints);
    std::vector<double> gradf2Values(numPoints);
    const double etaMax = this->kernelExtent();
    const double K1d = 0.5/simpsonsIntegration<volfunc<TableKernel<Dimension> >, double, double>(volfunc<TableKernel<Dimension> >(*this), 0.0, etaMax, numPoints);
    for (int i = 0; i < numPoints; ++i) {
      CHECK(i*mStepSize >= 0.0);
      const double zeta = i*mStepSize;
      f1Values[i] = f1Integral(*this, zeta, numPoints)/K1d;
      f2Values[i] = f2Integral(*this, zeta, numPoints)/K1d;
      // gradf1Values[i] = -f1Values[i]*f1Values[i]*gradf1Integral(*this, zeta, numPoints)*K1d;
    }
    // For now we do a kludgey numerical estimate of the gradient terms.
    for (int i = 0; i < numPoints; ++i) {
      const int i0 = max(i - 1, 0);
      const int i1 = min(i + 1, numPoints - 1);
      gradf1Values[i] = (f1Values[i] - f1Values[i0] + f1Values[i1] - f1Values[i])/((i1 - i0)*mStepSize);
      gradf2Values[i] = (f2Values[i] - f2Values[i0] + f2Values[i1] - f2Values[i])/((i1 - i0)*mStepSize);
    }
    setParabolicCoeffs(f1Values, mAf1, mBf1, mCf1);
    setParabolicCoeffs(f2Values, mAf2, mBf2, mCf2);
    setParabolicCoeffs(gradf1Values, mAgradf1, mBgradf1, mCgradf1);
    setParabolicCoeffs(gradf2Values, mAgradf2, mBgradf2, mCgradf2);
  }

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
    mAkernel = rhs.mAkernel;
    mBkernel = rhs.mBkernel;
    mCkernel = rhs.mCkernel;
    mAgrad = rhs.mAgrad;
    mBgrad = rhs.mBgrad;
    mCgrad = rhs.mCgrad;
    mAgrad2 = rhs.mAgrad2;
    mBgrad2 = rhs.mBgrad2;
    mCgrad2 = rhs.mCgrad2;
    mNumPoints = rhs.mNumPoints;
    mStepSize = rhs.mStepSize;
    mNperhValues = rhs.mNperhValues;
    mWsumValues =  rhs.mWsumValues;
    mMinNperh = rhs.mMinNperh;
    mMaxNperh = rhs.mMaxNperh;
    mAf1 = rhs.mAf1;
    mBf1 = rhs.mBf1;
    mCf1 = rhs.mCf1;
    mAf2 = rhs.mAf2;
    mBf2 = rhs.mBf2;
    mCf2 = rhs.mCf2;
  }
  return *this;
}

// //------------------------------------------------------------------------------
// // Linearly combine with another kernel.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// template<typename KernelType>
// void
// TableKernel<Dimension>::augment(const KernelType& kernel) {

//   // Set the volume normalization and kernel extent.
//   const double norm0 = this->volumeNormalization();
//   const double norm1 = kernel.volumeNormalization();
//   this->setVolumeNormalization(norm0 + norm1);
//   this->setInflectionPoint(0.5*(this->inflectionPoint() + kernel.inflectionPoint()));  // Punting for now.

//   // Fill in the kernel and gradient values.  Note that we will go ahead and fold
//   // the normalization constants in here, so we don't have to multiply by them
//   // in the value lookups.
//   const double deta = mStepSize;
//   for (int i = 0; i < mNumPoints; ++i) {
//     CHECK(i*mStepSize >= 0.0);
//     mKernelValues[i] = 0.5*(mKernelValues[i] + kernel(i*deta, 1.0));
//     mGradValues[i] = 0.5*(mGradValues[i] + kernel.grad(i*deta, 1.0));
//     mGrad2Values[i] = 0.5*(mGrad2Values[i] + kernel.grad2(i*deta, 1.0));
//   }

//   // Set the delta kernel values for internal use.
//   this->setParabolicCoeffs();

//   // // Adjust the table kernel values to reproduce the volume integral as closely as possible.
//   // const double vol1 = simpsonsVolumeIntegral<Dimension, TableKernel<Dimension> >(*this, 0.0, kernel.kernelExtent(), 10*numPoints);
//   // const double f = 1.0/vol1; // /(this->volumeNormalization());
//   // std::cerr << "Scaling: " << f << " " << vol1 << " " << this->volumeNormalization() << std::endl;
//   // for (int i = 0; i < numPoints; ++i) {
//   //   mKernelValues[i] *= f;
//   // }

//   // // Reset the delta kernel values for internal use.
//   // this->setDeltaKernelValues();

//   // Set the table of n per h values.
//   this->setNperhValues();

//   // That should be it, so we should have left the kernel in a valid state.
//   ENSURE(valid());
// }

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
setParabolicCoeffs(const std::vector<double>& table,
                   std::vector<double>& a,
                   std::vector<double>& b,
                   std::vector<double>& c) const {
  
  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> EMatrix;
  typedef Eigen::Matrix<double, 3, 1> EVector;

  // Because we do our fit in units of fixed eta step size, our A matrix of x elemements is constant and trivial!
  EMatrix A;
  A << 0.0, 0.0, 1.0,
       1.0, 1.0, 1.0,
       4.0, 2.0, 1.0;
  EMatrix Ainv = A.inverse();

  // Size stuff up.
  REQUIRE(mNumPoints > 0);
  a = std::vector<double>(mNumPoints);
  b = std::vector<double>(mNumPoints);
  c = std::vector<double>(mNumPoints);

  // Find the coefficient fits.
  for (int i0 = 0; i0 < mNumPoints - 2; ++i0) {
    const int i1 = i0 + 1;
    const int i2 = i0 + 2;
    CHECK(i2 < mNumPoints);
    EVector B, C;

    B << table[i0], table[i1], table[i2];
    C = Ainv*B;
    a[i1] = C(2);
    b[i1] = C(1);
    c[i1] = C(0);
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
  mWsumValues = vector<double>(mNumPoints);
  mNperhValues = vector<double>(mNumPoints);

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
  BEGIN_CONTRACT_SCOPE
  ENSURE(mWsumValues.size() == mNumPoints);
  ENSURE(mNperhValues.size() == mNumPoints);
  for (int i = 0; i != mNumPoints - 1; ++i) {
    ENSURE(mWsumValues[i] <= mWsumValues[i + 1]);
    ENSURE(mNperhValues[i] <= mNperhValues[i + 1]);
  }
  END_CONTRACT_SCOPE

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
          (int)mAkernel.size() == mNumPoints and
          (int)mBkernel.size() == mNumPoints and
          (int)mCkernel.size() == mNumPoints and
          (int)mAgrad.size() == mNumPoints and
          (int)mBgrad.size() == mNumPoints and
          (int)mCgrad.size() == mNumPoints and
          (int)mAgrad2.size() == mNumPoints and
          (int)mBgrad2.size() == mNumPoints and
          (int)mCgrad2.size() == mNumPoints and
          mStepSize > 0.0);
}

}
