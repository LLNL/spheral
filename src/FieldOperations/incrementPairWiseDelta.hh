#include "Utilities/rotationMatrix.hh"

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<typename Dimension>
void
incrementPairWiseDelta(typename Dimension::Vector& result,
                       typename Dimension::Vector& normalization,
                       const double& weight,
                       const double& fij,
                       const typename Dimension::Vector& rij,
                       const typename Dimension::SymTensor& Hi) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

  const Scalar frach = (0.01/Hi.Trace())/Dimension::nDim;
  const Vector rjiUnit = -rij.unitVector();
  const Tensor R = rotationMatrix(rjiUnit);
  const Tensor Rinverse = R.Transpose();
  const Scalar dxp = rij.magnitude()/(rij.magnitude2() + frach*frach);

  // Rotate the cumulative result to the prime frame.
  result = R*result;

  // Increment by the weighted pair wise contribtion in the
  // rotated frame.
  result.x(result.x() - weight*fij*dxp);

  // Rotate back and we're done.
  result = Rinverse*result;

  // Increment the normalization.
  typename Dimension::Vector dw;
  dw.x(weight);
  dw = Rinverse*dw;
  for (int i = 0; i < Dimension::nDim; ++i)
    normalization(i) += std::abs(dw(i));
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<typename Dimension>
void
incrementPairWiseDelta(typename Dimension::Tensor& result,
                       typename Dimension::Tensor& normalization,
                       const double& weight,
                       const typename Dimension::Vector& fij,
                       const typename Dimension::Vector& rij,
                       const typename Dimension::SymTensor& Hi) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

  const Scalar frach = (0.01/Hi.Trace())/Dimension::nDim;
  const Vector rjiUnit = -rij.unitVector();
  const Tensor R = rotationMatrix(rjiUnit);
  const Tensor Rinverse = R.Transpose();
  const Scalar dxp = rij.magnitude()/(rij.magnitude2() + frach*frach);

  // Represent this pairs delta as the first column of a matrix
  // in the rotated frame.
  const Vector dfdxpcol = -dxp*R*fij;
  Tensor dfdxp;
  dfdxp.setColumn(0, dfdxpcol);

  // Rotate the delta to the lab frame.
  dfdxp.rotationalTransform(Rinverse);

  // Increment the result with the delta.
  result += weight*dfdxp;

  // Update the normalization.
  Tensor tweight(0.0);
  for (int i = 0; i < Dimension::nDim; ++i) tweight(i,0) = weight;
  tweight.rotationalTransform(Rinverse);
  for (int i = 0; i < Dimension::nDim; ++i) {
    for (int j = 0; j < Dimension::nDim; ++j) {
      normalization(i,j) += std::abs(tweight(i,j));
    }
  }
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<typename Dimension>
void
renormalizeGradient(typename Dimension::Vector& result,
                    const typename Dimension::Vector& normalization) {
  for (int i = 0; i < Dimension::nDim; ++i) {
    result(i) /= normalization(i) + 1.0e-10;
  }
}

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<typename Dimension>
void
renormalizeGradient(typename Dimension::Tensor& result,
                    const typename Dimension::Tensor& normalization) {
  for (int i = 0; i < Dimension::nDim; ++i) {
    for (int j = 0; j < Dimension::nDim; ++j) {
      result(i,j) /= normalization(i,j) + 1.0e-10;
    }
  }
}

