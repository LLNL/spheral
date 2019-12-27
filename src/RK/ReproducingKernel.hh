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

  // Constructor
  explicit ReproducingKernel(const TableKernel<Dimension>& W,
                             const RKOrder order);
                             
  ~ReproducingKernel();

  // Base kernel calls
  Scalar evaluateBaseKernel(const Vector& x,
                            const SymTensor& H);
  Vector evaluateBaseGradient(const Vector& x,
                              const SymTensor& H);
  SymTensor evaluateBaseHessian(const Vector& x,
                                const SymTensor& H);
  std::pair<Scalar, Vector> evaluateBaseKernelAndGradient(const Vector& x,
                                                          const SymTensor& H);

  // Corrected kernel calls
  Scalar evaluateKernel(const Vector& x,
                        const SymTensor& H,
                        const std::vector<double>& corrections);
  Vector evaluateGradient(const Vector& x,
                          const SymTensor& H,
                          const std::vector<double>& corrections);
  SymTensor evaluateHessian(const Vector& x,
                            const SymTensor& H,
                            const std::vector<double>& corrections);
  std::pair<Scalar, Vector> evaluateKernelAndGradient(const Vector& x,
                                                      const SymTensor& H,
                                                      const std::vector<double>& corrections);
  std::tuple<Scalar, Vector, Scalar> evaluateKernelAndGradients(const Vector& x,
                                                                const SymTensor& H,
                                                                const std::vector<double>& corrections);

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mW;   // The base interpolation kernel
  RKOrder mOrder;

  // Pointers to the correct member methods of RKUtilities
  Scalar                    (*mEvaluateBaseKernel)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  Scalar                    (*mEvaluateBaseGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  Scalar                    (*mEvaluateBaseHessian)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  std::pair<Scalar, Vector> (*mEvaluateBaseKernelAndGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);

  Scalar                             (*mEvaluateKernel)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  Scalar                             (*mEvaluateGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  Scalar                             (*mEvaluateHessian)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  std::pair<Scalar, Vector>          (*mEvaluateKernelAndGradient)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);
  std::tuple<Scalar, Vector, Scalar> (*mEvaluateKernelAndGradients)(const TableKernel<Dimension>&, const Vector&, const SymTensor&);

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

  // No default constructor
  ReproducingKernel();
};

}

#include "RK/ReproducingKernelInline.hh"

#endif
