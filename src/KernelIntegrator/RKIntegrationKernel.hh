//---------------------------------Spheral++----------------------------------//
// RKIntegrationKernel
//----------------------------------------------------------------------------//
#ifndef __Spheral_RKIntegrationKernel__
#define __Spheral_RKIntegrationKernel__

#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"
#include "IntegrationKernel.hh"
#include "SPHIntegrationKernel.hh"

#include <utility>
#include <vector>

namespace Spheral {

template<typename Dimension, int order>
class RKIntegrationKernel : public IntegrationKernel<Dimension> {
public:
  // Static size information for arrays
  static constexpr int getPolynomialSize(int d, int o) {
    return (d == 1 ? o + 1
            : d == 2 ? (1 + o) * (2 + o) / 2
            : (1 + o) * (2 + o) * (3 + o) / 6);
  }
  static constexpr int polynomialSize = getPolynomialSize(Dimension::nDim, order);
  static constexpr int gradPolynomialSize = polynomialSize * Dimension::nDim;
  static constexpr int correctionsSize = polynomialSize * (1 + Dimension::nDim);
  static constexpr std::array<int, Dim<3>::nDim> offsetGradC = {polynomialSize, 2 * polynomialSize, 3 * polynomialSize};
  static constexpr std::array<int, Dim<3>::nDim> offsetGradP = {0, polynomialSize, 2 * polynomialSize};
  
  // Typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename std::array<Scalar, polynomialSize> PolyArray;
  typedef typename std::array<Scalar, gradPolynomialSize> GradPolyArray;
  typedef Eigen::Matrix<Scalar, polynomialSize, 1> EigenVector;
  typedef Eigen::Matrix<Scalar, polynomialSize, polynomialSize> EigenMatrix;
  typedef std::array<EigenVector, Dimension::nDim> VectorOfEigenVector;
  typedef std::array<EigenMatrix, Dimension::nDim> VectorOfEigenMatrix;
  typedef std::array<std::array<Scalar, Dimension::nDim>, order + 1> PolyArray1d;
  typedef std::vector<Scalar> CorrectionsVector;
  
  RKIntegrationKernel(const TableKernel<Dimension>& kernel);
  
  virtual double extent(const Scalar Hmult) const { return mKernel.kernelExtent() / Hmult; }
  
  // Evaluate the all the functions at the point xp
  virtual void evaluate(const Vector& xp,
                        const std::vector<std::pair<int, int>>& indices,
                        const FieldList<Dimension, Vector>& position,
                        const FieldList<Dimension, SymTensor>& H,
                        const FieldList<Dimension, Scalar>& volume,
                        const Scalar Hmult,
                        std::vector<Scalar>& values,
                        std::vector<Vector>& dvalues) const override;
  
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Need this until C++17
  
protected:

  // Get the polynomial vectors
  void getPolynomials(const Vector& x,
                      PolyArray& p,
                      GradPolyArray& dp) const;
  
  // Helper function that calculates the corrections
  void corrections(const Vector& xp,
                   const std::vector<std::pair<int, int>>& indices,
                   const FieldList<Dimension, Vector>& position,
                   const FieldList<Dimension, Scalar>& volume,
                   const std::vector<Scalar>& values,
                   const std::vector<Vector>& dvalues,
                   CorrectionsVector& corrections) const;

  // Helper function that replaces the SPH values by RK values
  void replace(const Vector& xp,
               const std::vector<std::pair<int, int>>& indices,
               const FieldList<Dimension, Vector>& position,
               const CorrectionsVector& corrections,
               std::vector<Scalar>& values,
               std::vector<Vector>& dvalues) const;
  
  // Do inner product, given offsets
  template<typename DataType1, typename DataType2>
  static inline Scalar innerProductRK(const DataType1& x,
                                      const DataType2& y,
                                      const int offsetx,
                                      const int offsety);
  
  // Input data
  const TableKernel<Dimension>& mKernel;

  // Standard SPH kernel
  SPHIntegrationKernel<Dimension> mSPHKernel;
  
  // Scratch variables
  mutable EigenMatrix mM;          // poly*poly matrices
  mutable VectorOfEigenMatrix mDM;  
  mutable EigenVector mC;          // corrections vectors
  mutable VectorOfEigenVector mDC;
  mutable EigenVector mRHS;        // for solving matrix equations
  mutable PolyArray mP;            // poly vectors
  mutable GradPolyArray mDP;
  mutable CorrectionsVector mCV; // combined corrections vector
  mutable PolyArray1d mQ;          // for evaluating polynomials of high order
  mutable PolyArray1d mDQ;
};

} // end namespace Spheral

#include "RKIntegrationKernelInline.hh"

#endif
