//---------------------------------Spheral++----------------------------------//
// ReproducingKernel
//
// Front-end class to wrap the Reproducing Kernel corrected methods
//----------------------------------------------------------------------------//
#ifndef __LLNLSpheral_ReproducingKernel__
#define __LLNLSpheral_ReproducingKernel__

#include "RK/RKUtilities.hh"

namespace Spheral {

template<typename Dimension>
class ReproducingKernel {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructor
  ReproducingKernel(const TableKernel<Dimension>& W,
                    const RKOrder order);
  ReproducingKernel();
  ReproducingKernel(const ReproducingKernel& rhs);
  ReproducingKernel& operator=(const ReproducingKernel& rhs);
  ~ReproducingKernel();

  // Base kernel calls
  Scalar evaluateBaseKernel(const Vector& x,
                            const SymTensor& H) const;
  Vector evaluateBaseGradient(const Vector& x,
                              const SymTensor& H) const;
  SymTensor evaluateBaseHessian(const Vector& x,
                                const SymTensor& H) const;
  std::pair<Scalar, Vector> evaluateBaseKernelAndGradient(const Vector& x,
                                                          const SymTensor& H) const;

  // Corrected kernel calls
  Scalar evaluateKernel(const Vector& x,
                        const SymTensor& H,
                        const std::vector<double>& corrections) const;
  Vector evaluateGradient(const Vector& x,
                          const SymTensor& H,
                          const std::vector<double>& corrections) const;
  SymTensor evaluateHessian(const Vector& x,
                            const SymTensor& H,
                            const std::vector<double>& corrections) const;
  std::pair<Scalar, Vector> evaluateKernelAndGradient(const Vector& x,
                                                      const SymTensor& H,
                                                      const std::vector<double>& corrections) const;
  std::tuple<Scalar, Vector, Scalar> evaluateKernelAndGradients(const Vector& x,
                                                                const SymTensor& H,
                                                                const std::vector<double>& corrections) const;

  // Compute corrections and normals
  void computeCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                          const FieldList<Dimension, Scalar>& volume,
                          const FieldList<Dimension, Vector>& position,
                          const FieldList<Dimension, SymTensor>& H,
                          const bool needHessian,
                          FieldList<Dimension, std::vector<double>>& zerothCorrections,
                          FieldList<Dimension, std::vector<double>>& corrections) const;
  void computeNormal(const ConnectivityMap<Dimension>& connectivityMap,
                     const FieldList<Dimension, Scalar>& volume,
                     const FieldList<Dimension, Vector>& position,
                     const FieldList<Dimension, SymTensor>& H,
                     const FieldList<Dimension, std::vector<double>>& corrections,
                     FieldList<Dimension, Scalar>& surfaceArea,
                     FieldList<Dimension, Vector>& normal) const;

  // Apply a transformation operator to a corrections vector
  void applyTransformation(const typename Dimension::Tensor& T,
                           std::vector<double>& corrections) const;

  // Access the internal state
  RKOrder order() const;
  const TableKernel<Dimension>& kernel() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>* mWptr;   // The base interpolation kernel
  RKOrder mOrder;

  // Pointers to the correct member methods of RKUtilities
  Scalar                    (*mEvaluateBaseKernel)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  Vector                    (*mEvaluateBaseGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  SymTensor                 (*mEvaluateBaseHessian)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  std::pair<Scalar, Vector> (*mEvaluateBaseKernelAndGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);

  Scalar                             (*mEvaluateKernel)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const std::vector<double>&);
  Vector                             (*mEvaluateGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const std::vector<double>&);
  SymTensor                          (*mEvaluateHessian)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const std::vector<double>&);
  std::pair<Scalar, Vector>          (*mEvaluateKernelAndGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const std::vector<double>&);
  std::tuple<Scalar, Vector, Scalar> (*mEvaluateKernelAndGradients)(const TableKernel<Dimension>&, const Vector&, const SymTensor&, const std::vector<double>&);

  void (*mComputeCorrections)(const ConnectivityMap<Dimension>&,
                              const TableKernel<Dimension>&,
                              const FieldList<Dimension, Scalar>&,
                              const FieldList<Dimension, Vector>&,
                              const FieldList<Dimension, SymTensor>&,
                              const bool,
                              FieldList<Dimension, std::vector<double>>&,
                              FieldList<Dimension, std::vector<double>>&);
  void (*mComputeNormal)(const ConnectivityMap<Dimension>&,
                         const TableKernel<Dimension>&,
                         const FieldList<Dimension, Scalar>&,
                         const FieldList<Dimension, Vector>&,
                         const FieldList<Dimension, SymTensor>&,
                         const FieldList<Dimension, std::vector<double>>&,
                         FieldList<Dimension, Scalar>&,
                         FieldList<Dimension, Vector>&);
  void (*mApplyTransformation)(const typename Dimension::Tensor& T,
                               std::vector<double>& corrections);
};

}

#include "RK/ReproducingKernelInline.hh"

#endif
