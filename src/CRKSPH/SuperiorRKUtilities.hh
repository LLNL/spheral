//---------------------------------Spheral++----------------------------------//
// SuperiorRKUtilities
//
// Computes and evaluates RK corrections
//----------------------------------------------------------------------------//
#ifndef __LLNLSpheral_SuperiorRKCorrections__
#define __LLNLSpheral_SuperiorRKCorrections__

#include <vector>
#include "CRKSPH/CRKSPHCorrectionParams.hh"
#include "Field/FieldList.hh"

template<typename Dimension, CRKOrder correctionOrder>
class SuperiorRKUtilities {
public:
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

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
  template<typename DataType> static FieldList<Dimension, typename DataType>
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
  
  // Get starting index for coefficient array derivatives
  static inline int offsetGradC(const int d);
  static inline int offsetHessC(const int d1, const int d2);

  // Get starting index for gradPolynomial and hessPolynomial
  static inline int offsetGradP(const int d);
  static inline int offsetHessP(const int d1, const int d2);
  
  // Get the polynomial vectors
  static inline std::vector<double> getPolynomials(Vector& x);
  static inline std::vector<double> getGradPolynomials(Vector& x);
  static inline std::vector<double> getHessPolynomials(Vector& x);

  // The size of the polynomials and coefficients, not including derivatives
  static int polynomialSize;
};

#include "SuperiorRKUtilitiesInline.hh"
  
#endif
