//---------------------------------Spheral++----------------------------------//
// RKIntegrationKernel
//----------------------------------------------------------------------------//
#include "RKIntegrationKernel.hh"

#include "Eigen/Dense"

namespace Spheral {

//------------------------------------------------------------------------------
// We can remove these as of C++17
//------------------------------------------------------------------------------
template<typename Dimension, int order>
constexpr std::array<int, Dim<3>::nDim>
RKIntegrationKernel<Dimension, order>::offsetGradC;

template<typename Dimension, int order>
constexpr std::array<int, Dim<3>::nDim>
RKIntegrationKernel<Dimension, order>::offsetGradP;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension, int order>
RKIntegrationKernel<Dimension, order>::
RKIntegrationKernel(const TableKernel<Dimension>& kernel):
  mKernel(kernel),
  mSPHKernel(kernel),
  mCV(correctionsSize) {
}

//------------------------------------------------------------------------------
// Evaluate all the functions at a point
//------------------------------------------------------------------------------
template<typename Dimension, int order>
void
RKIntegrationKernel<Dimension, order>::
evaluate(const Vector& xp,
         const std::vector<std::pair<int, int>>& indices,
         const FieldList<Dimension, Vector>& position,
         const FieldList<Dimension, SymTensor>& H,
         const FieldList<Dimension, Scalar>& volume,
         const Scalar Hmult,
         std::vector<Scalar>& values,
         std::vector<Vector>& dvalues) const {
  CHECK(values.size() == indices.size());
  CHECK(dvalues.size() == indices.size());
  
  // Get the base values for the kernel
  mSPHKernel.evaluate(xp, indices, position, H, volume, Hmult,
                      values, dvalues);

  // Get the corrections
  corrections(xp, indices, position, volume, values, dvalues,
              mCV);

  // Evaluate the values using the corrections
  replace(xp, indices, position, mCV,
          values, dvalues);
}

//------------------------------------------------------------------------------
// Compute the corrections for a point, given the values
//------------------------------------------------------------------------------
template<typename Dimension, int order>
void
RKIntegrationKernel<Dimension, order>::
corrections(const Vector& xp,
            const std::vector<std::pair<int, int>>& indices,
            const FieldList<Dimension, Vector>& position,
            const FieldList<Dimension, Scalar>& volume,
            const std::vector<Scalar>& values,
            const std::vector<Vector>& dvalues,
            CorrectionsVector& corrections) const {
  // Initialize the matrices
  mM.setZero();
  for (auto d = 0; d < Dimension::nDim; ++d) {
    mDM[d].setZero();
  }
  
  // Sum up the matrix quantities
  const auto numIndices = indices.size();
  CHECK(values.size() == numIndices);
  CHECK(dvalues.size() == numIndices);
  for (auto j = 0u; j < numIndices; ++j) {
    const auto nodeListj = indices[j].first;
    const auto nodej = indices[j].second;
    const auto xj = position(nodeListj, nodej);
    const auto vj = volume(nodeListj, nodej);
    const auto xpj = xp - xj;
    getPolynomials(xpj, mP, mDP);
    for (auto k = 0; k < polynomialSize; ++k) {
      for (auto l = k; l < polynomialSize; ++l) {
        mM(k, l) += vj * mP[k] * mP[l] * values[j];
        
        for (auto d = 0; d < Dimension::nDim; ++d) {
          const auto offP = offsetGradP[d];
          mDM[d](k, l) += vj * ((mDP[offP+k] * mP[l] + mP[k] * mDP[offP+l]) * values[j] + mP[k] * mP[l] * dvalues[j](d));
        }
      }
    }
  }
  
  // Fill in the matrix symmetries
  for (auto k = 0; k < polynomialSize; ++k) {
    for (auto l = 0; l < k; ++l) {
      mM(k, l) = mM(l, k);
      for (auto d = 0; d < Dimension::nDim; ++d) {
        mDM[d](k, l) = mDM[d](l, k);
      }
    }
  }

  // Get inverse of the M matrix
  auto solver = mM.colPivHouseholderQr();
  
  // Calculate the corrections
  mRHS.setZero();
  mRHS(0) = 1;
  mC = solver.solve(mRHS);

  // Calculate the gradient corrections
  for (auto d = 0; d < Dimension::nDim; ++d) {
    mRHS = -(mDM[d] * mC);
    mDC[d] = solver.solve(mRHS);
  }

  // We're going to keep the corrections a vector for now, which is a bit of
  // extra effort here but keeps compatibility with RKUtilities
  corrections.resize(correctionsSize);
  for (auto k = 0; k < polynomialSize; ++k) {
    corrections[k] = mC(k);
    for (auto d = 0; d < Dimension::nDim; ++d) {
      const auto offd = offsetGradC[d];
      corrections[offd+k] = mDC[d](k);
    }
  }
}

//------------------------------------------------------------------------------
// Replace the SPH values by RK values
//------------------------------------------------------------------------------
template<typename Dimension, int order>
void
RKIntegrationKernel<Dimension, order>::
replace(const Vector& xp,
        const std::vector<std::pair<int, int>>& indices,
        const FieldList<Dimension, Vector>& position,
        const CorrectionsVector& corrections,
        std::vector<Scalar>& values,
        std::vector<Vector>& dvalues) const {
  const auto numIndices = indices.size();
  CHECK(values.size() == numIndices);
  CHECK(dvalues.size() == numIndices);
  for (auto j = 0u; j < numIndices; ++j) {
    const auto nodeListj = indices[j].first;
    const auto nodej = indices[j].second;
    const auto xj = position(nodeListj, nodej);
    const auto xpj = xp - xj;
    
    getPolynomials(xpj, mP, mDP);
    
    // We need to evaluate the derivative first, since it depends on the kernel
    const auto CP = innerProductRK(corrections, mP, 0, 0);
    for (auto d = 0; d < Dimension::nDim; ++d) {
      const auto CdP = innerProductRK(corrections, mDP, 0, offsetGradP[d]);
      const auto dCP = innerProductRK(corrections, mP, offsetGradC[d], 0);
      dvalues[j](d) = (CdP + dCP) * values[j] + CP * dvalues[j](d);
    }
    
    // Now we can overwrite the kernel
    values[j] = CP * values[j];
  }
}

} // end namespace Spheral
