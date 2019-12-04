//---------------------------------Spheral++----------------------------------//
// RKUtilities
//
// Computes and evaluates RK corrections
//----------------------------------------------------------------------------//
#ifndef __LLNLSpheral_RKUtilities__
#define __LLNLSpheral_RKUtilities__

#include <vector>
#include "CRKSPH/CRKSPHCorrectionParams.hh"
#include "Field/FieldList.hh"

namespace Spheral {

template<typename Dimension, CRKOrder correctionOrder>
class RKUtilities {
public:
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Get the polynomial vectors
  static inline std::vector<double> getPolynomials(const Vector& x);
  static inline std::vector<double> getGradPolynomials(const Vector& x);
  static inline std::vector<double> getHessPolynomials(const Vector& x);
  
  // Evaluate base functions
  static Scalar evaluateBaseKernel(const TableKernel<Dimension>& kernel,
                                   const Vector& x,
                                   const SymTensor& H);
  static Vector evaluateBaseGradient(const TableKernel<Dimension>& kernel,
                                     const Vector& x,
                                     const SymTensor& H);
  static SymTensor evaluateBaseHessian(const TableKernel<Dimension>& kernel,
                                       const Vector& x,
                                       const SymTensor& H);

  // Evaluate kernels
  static Scalar evaluateKernel(const TableKernel<Dimension>& kernel,
                               const Vector& x,
                               const SymTensor& H,
                               const std::vector<double>& corrections);
  static Vector evaluateGradient(const TableKernel<Dimension>& kernel,
                                 const Vector& x,
                                 const SymTensor& H,
                                 const std::vector<double>& corrections);
  static SymTensor evaluateHessian(const TableKernel<Dimension>& kernel,
                                   const Vector& x,
                                   const SymTensor& H,
                                   const std::vector<double>& corrections);
  
  // Compute the corrections
  static void computeCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                                 const TableKernel<Dimension>& kernel,
                                 const FieldList<Dimension, Scalar>& volume,
                                 const FieldList<Dimension, Vector>& position,
                                 const FieldList<Dimension, SymTensor>& H,
                                 const bool needHessian,
                                 FieldList<Dimension, std::vector<double>>& corrections);

  // Interpolate a field
  template<typename DataType> static FieldList<Dimension, DataType>
  interpolateField(const TableKernel<Dimension>& kernel,
                   const FieldList<Dimension, Scalar>& volume,
                   const FieldList<Dimension, Vector>& position,
                   const FieldList<Dimension, SymTensor>& H,
                   const FieldList<Dimension, std::vector<double>>& corrections,
                   const bool needHessian,
                   const FieldList<Dimension, DataType>& field,
                   FieldList<Dimension, DataType>& interpolant);
  template<typename DataType> static FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
  gradientField(const TableKernel<Dimension>& kernel,
                const FieldList<Dimension, Scalar>& volume,
                const FieldList<Dimension, Vector>& position,
                const FieldList<Dimension, SymTensor>& H,
                const FieldList<Dimension, std::vector<double>>& corrections,
                const bool needHessian,
                const FieldList<Dimension, DataType>& field,
                FieldList<Dimension, DataType>& interpolant);

  // Do inner product, given offsets
  static inline Scalar innerProductRK(const std::vector<double>& x,
                                      const std::vector<double>& y,
                                      const int offsetx,
                                      const int offsety);

  // Get expected size of corrections
  static inline int correctionsSize(const bool needHessian);
  
  // Get storage size of a symmetric matrix
  static inline int symmetricMatrixSize(const int d);
  
  // Get flat index for a symmetric set of indices
  static inline int flatSymmetricIndex(const int d1, const int d2);
  
  // Get starting index for coefficient array derivatives
  static inline int offsetGradC(const int d);
  static inline int offsetHessC(const int d1, const int d2);

  // Get starting index for gradPolynomial and hessPolynomial
  static inline int offsetGradP(const int d);
  static inline int offsetHessP(const int d1, const int d2);
  
  // The size of the polynomials and coefficients, not including derivatives
  static constexpr int polynomialSize = (correctionOrder == CRKOrder::ZerothOrder ? 1
                                         : (correctionOrder == CRKOrder::LinearOrder
                                            ? (Dimension::nDim == 1 ? 2 : (Dimension::nDim == 2 ? 3 : 4))
                                            : (correctionOrder == CRKOrder::QuadraticOrder
                                               ? (Dimension::nDim == 1 ? 3 : (Dimension::nDim == 2 ? 6 : 10))
                                               : (correctionOrder == CRKOrder::CubicOrder
                                                  ? (Dimension::nDim == 1 ? 4 : (Dimension::nDim == 2 ? 10 : 20))
                                                  : (correctionOrder == CRKOrder::QuarticOrder
                                                     ? (Dimension::nDim == 1 ? 5 : (Dimension::nDim == 2 ? 15 : 35))
                                                     : (correctionOrder == CRKOrder::QuinticOrder
                                                        ? (Dimension::nDim == 1 ? 6 : (Dimension::nDim == 2 ? 21 : 56))
                                                        : (correctionOrder == CRKOrder::SexticOrder
                                                           ? (Dimension::nDim == 1 ? 7 : (Dimension::nDim == 2 ? 28 : 84))
                                                           : (correctionOrder == CRKOrder::SepticOrder
                                                              ? (Dimension::nDim == 1 ? 8 : (Dimension::nDim == 2 ? 36 : 120))
                                                              : -1)))))))); // if order not found, return -1 to produce error
};

} // end namespace Spheral

#include "RKUtilitiesInline.hh"
  
#endif
