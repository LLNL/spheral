//---------------------------------Spheral++----------------------------------//
// RKUtilities
//
// Computes and evaluates RK corrections
//----------------------------------------------------------------------------//
#ifndef __LLNLSpheral_RKUtilities__
#define __LLNLSpheral_RKUtilities__

#include "RKCorrectionParams.hh"
#include "Field/FieldList.hh"
#include <vector>
#include <tuple>

namespace Spheral {

template<typename Dimension, RKOrder correctionOrder>
class RKUtilities {
public:
  // Get storage size of a symmetric matrix
  static constexpr int symmetricMatrixSize(int d) {
    return d * (d + 1) / 2;
  }
  
  // // Size of the sparse storage for the Hessian
  // static constexpr int hessBaseSize = symmetricMatrixSize(Dimension::nDim);//(Dimension::nDim == 1 ? 1 : (Dimension::nDim == 2 ? 3 : 6));
  
  // The size of the polynomial arrays
  static constexpr int polynomialSize = (correctionOrder == RKOrder::ZerothOrder ? 1 
                                         : correctionOrder == RKOrder::LinearOrder
                                         ? (Dimension::nDim == 1 ? 2 : (Dimension::nDim == 2 ? 3 : 4))
                                         : correctionOrder == RKOrder::QuadraticOrder
                                         ? (Dimension::nDim == 1 ? 3 : (Dimension::nDim == 2 ? 6 : 10))
                                         : correctionOrder == RKOrder::CubicOrder
                                         ? (Dimension::nDim == 1 ? 4 : (Dimension::nDim == 2 ? 10 : 20))
                                         : correctionOrder == RKOrder::QuarticOrder
                                         ? (Dimension::nDim == 1 ? 5 : (Dimension::nDim == 2 ? 15 : 35))
                                         : correctionOrder == RKOrder::QuinticOrder
                                         ? (Dimension::nDim == 1 ? 6 : (Dimension::nDim == 2 ? 21 : 56))
                                         : correctionOrder == RKOrder::SexticOrder
                                         ? (Dimension::nDim == 1 ? 7 : (Dimension::nDim == 2 ? 28 : 84))
                                         : correctionOrder == RKOrder::SepticOrder
                                         ? (Dimension::nDim == 1 ? 8 : (Dimension::nDim == 2 ? 36 : 120))
                                         : -1); // if order not found, return -1 to produce error
  static constexpr int gradPolynomialSize = polynomialSize * Dimension::nDim;
  static constexpr int hessPolynomialSize = polynomialSize * symmetricMatrixSize(Dimension::nDim);

  // Due to the hessian being optional, the size of the corrections vector in a dimension is not unique
  // The below integer is added to the hessian size for certain orders to make the hessian size unique
  // The duplicates occur between o2grad/o1hess for all dimensions, and o5grad/o3hess in 1D
  static constexpr int hessUniqueSizeModification = (correctionOrder == RKOrder::LinearOrder
                                                     || (Dimension::nDim == 1 && correctionOrder == RKOrder::CubicOrder)
                                                     ? 1
                                                     : 0);
  
  // The size of the corrections, with kernel + gradient only and with kernel + gradient + hessian
  static constexpr int gradCorrectionsSize = polynomialSize * (1 + Dimension::nDim);
  static constexpr int hessCorrectionsSize = polynomialSize * (1 + Dimension::nDim + symmetricMatrixSize(Dimension::nDim)) + hessUniqueSizeModification;
  
  // Typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename std::array<double, polynomialSize> PolyArray;
  typedef typename std::array<double, gradPolynomialSize> GradPolyArray;
  typedef typename std::array<double, hessPolynomialSize> HessPolyArray;

  // Get the polynomial vectors
  static inline void getPolynomials(const Vector& x,
                                    PolyArray& p);
  static inline void getGradPolynomials(const Vector& x,
                                        GradPolyArray& p);
  static inline void getHessPolynomials(const Vector& x,
                                        HessPolyArray& p);
  
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
  static std::pair<Scalar, Vector> evaluateBaseKernelAndGradient(const TableKernel<Dimension>& kernel,
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
  static std::pair<Scalar, Vector> evaluateKernelAndGradient(const TableKernel<Dimension>& kernel,
                                                             const Vector& x,
                                                             const SymTensor& H,
                                                             const std::vector<double>& corrections);

  // This one returns the (RK kernel, RK gradient, base gradient magnitude)
  static std::tuple<Scalar, Vector, Scalar> evaluateKernelAndGradients(const TableKernel<Dimension>& kernel,
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
                                 FieldList<Dimension, std::vector<double>>& zerothCorrections,
                                 FieldList<Dimension, std::vector<double>>& corrections);

  // Get a guess for the surface normals - best if done with zeroth order
  static void computeNormal(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& kernel,
                            const FieldList<Dimension, Scalar>& volume,
                            const FieldList<Dimension, Vector>& position,
                            const FieldList<Dimension, SymTensor>& H,
                            const FieldList<Dimension, std::vector<double>>& corrections,
                            FieldList<Dimension, Scalar>& surfaceArea,
                            FieldList<Dimension, Vector>& normal);
  
  // // Interpolate a field
  // template<typename DataType> static FieldList<Dimension, DataType>
  // interpolateField(const TableKernel<Dimension>& kernel,
  //                  const FieldList<Dimension, Scalar>& volume,
  //                  const FieldList<Dimension, Vector>& position,
  //                  const FieldList<Dimension, SymTensor>& H,
  //                  const FieldList<Dimension, std::vector<double>>& corrections,
  //                  const bool needHessian,
  //                  const FieldList<Dimension, DataType>& field,
  //                  FieldList<Dimension, DataType>& interpolant);
  // template<typename DataType> static FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
  // gradientField(const TableKernel<Dimension>& kernel,
  //               const FieldList<Dimension, Scalar>& volume,
  //               const FieldList<Dimension, Vector>& position,
  //               const FieldList<Dimension, SymTensor>& H,
  //               const FieldList<Dimension, std::vector<double>>& corrections,
  //               const bool needHessian,
  //               const FieldList<Dimension, DataType>& field,
  //               FieldList<Dimension, DataType>& interpolant);

  // Do inner product, given offsets
  template<typename DataType1, typename DataType2>
  static inline Scalar innerProductRK(const DataType1& x,
                                      const DataType2& y,
                                      const int offsetx,
                                      const int offsety);

  // Get expected size of corrections
  static inline int correctionsSize(const bool needHessian); 
  static inline int zerothCorrectionsSize(const bool needHessian);
  
  // Get flat index for a symmetric set of indices
  static inline int flatSymmetricIndex(const int d1, const int d2);
  
  // Get starting index for coefficient array derivatives
  static inline int offsetGradC(const int d);
  static inline int offsetHessC(const int d1, const int d2);

  // Get starting index for gradPolynomial and hessPolynomial
  static inline int offsetGradP(const int d);
  static inline int offsetHessP(const int d1, const int d2);
};

//------------------------------------------------------------------------------
// Provide frontends to workaround templating of correction order
//------------------------------------------------------------------------------
// Get correction size, given the correction order and whether the hessian is needed
template<typename Dimension>
int
RKCorrectionsSize(const RKOrder correctionOrder,
                  const bool needHessian);

// Get correction order and whether the hessian is used, given the corrections size
template<typename Dimension>
std::pair<RKOrder, bool>
RKOrderAndHessian(const int size);

// RK kernel
template<typename Dimension>
typename Dimension::Scalar
RKKernel(const TableKernel<Dimension>& W,
         const typename Dimension::Vector& x,
         const typename Dimension::SymTensor& H,
         const std::vector<double>& corrections);
         
// RK gradient
template<typename Dimension>
typename Dimension::Vector
RKGradient(const TableKernel<Dimension>& W,
           const typename Dimension::Vector& x,
           const typename Dimension::SymTensor& H,
           const std::vector<double>& corrections);

// RK Hessian
template<typename Dimension>
typename Dimension::SymTensor
RKHessian(const TableKernel<Dimension>& W,
          const typename Dimension::Vector& x,
          const typename Dimension::SymTensor& H,
          const std::vector<double>& corrections);

// RK kernel + gradient
template<typename Dimension>
void
RKKernelAndGradient(typename Dimension::Scalar& WRK,
                    typename Dimension::Vector& gradWSPH,
                    typename Dimension::Vector& gradWRK,
                    const TableKernel<Dimension>& W,
                    const typename Dimension::Vector& x,
                    const typename Dimension::SymTensor& H,
                    const std::vector<double>& corrections);

// RK base kernel
template<typename Dimension>
typename Dimension::Scalar
RKBaseKernel(const TableKernel<Dimension>& W,
             const typename Dimension::Vector& x,
             const typename Dimension::SymTensor& H);
         
// RK base gradient
template<typename Dimension>
typename Dimension::Vector
RKBaseGradient(const TableKernel<Dimension>& W,
               const typename Dimension::Vector& x,
               const typename Dimension::SymTensor& H);

// RK base Hessian
template<typename Dimension>
typename Dimension::SymTensor
RKBaseHessian(const TableKernel<Dimension>& W,
              const typename Dimension::Vector& x,
              const typename Dimension::SymTensor& H);

// Compute corrections
template<typename Dimension>
void
computeRKCorrections(const RKOrder order,
                     const ConnectivityMap<Dimension>& connectivityMap,
                     const TableKernel<Dimension>& kernel,
                     const FieldList<Dimension, typename Dimension::Scalar>& volume,
                     const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     const bool needHessian,
                     FieldList<Dimension, std::vector<double>>& zerothCorrections,
                     FieldList<Dimension, std::vector<double>>& corrections);

// Surface normals
template<typename Dimension>
void
computeRKNormal(const RKOrder order,
                const ConnectivityMap<Dimension>& connectivityMap,
                const TableKernel<Dimension>& kernel,
                const FieldList<Dimension, typename Dimension::Scalar>& volume,
                const FieldList<Dimension, typename Dimension::Vector>& position,
                const FieldList<Dimension, typename Dimension::SymTensor>& H,
                const FieldList<Dimension, std::vector<double>>& corrections,
                FieldList<Dimension, typename Dimension::Scalar>& surfaceArea,
                FieldList<Dimension, typename Dimension::Vector>& normal);
} // end namespace Spheral

#include "RKUtilitiesInline.hh"
  
#endif
