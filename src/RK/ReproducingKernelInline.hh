namespace Spheral {

//------------------------------------------------------------------------------
// evaluateBaseKernel
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ReproducingKernel<Dimension>::
evaluateBaseKernel(const typename Dimension::Vector& x,
                   const typename Dimension::SymTensor& H) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateBaseKernel)(*mWptr, x, H);
}

//------------------------------------------------------------------------------
// evaluateBaseGradient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
ReproducingKernel<Dimension>::
evaluateBaseGradient(const typename Dimension::Vector& x,
                     const typename Dimension::SymTensor& H) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateBaseGradient)(*mWptr, x, H);
}

//------------------------------------------------------------------------------
// evaluateBaseHessian
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::SymTensor
ReproducingKernel<Dimension>::
evaluateBaseHessian(const typename Dimension::Vector& x,
                    const typename Dimension::SymTensor& H) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateBaseHessian)(*mWptr, x, H);
}

//------------------------------------------------------------------------------
// evaluateBaseKernelAndGradient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::pair<typename Dimension::Scalar, typename Dimension::Vector>
ReproducingKernel<Dimension>::
evaluateBaseKernelAndGradient(const typename Dimension::Vector& x,
                     const typename Dimension::SymTensor& H) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateBaseKernelAndGradient)(*mWptr, x, H);
}

//------------------------------------------------------------------------------
// evaluateKernel
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
ReproducingKernel<Dimension>::
evaluateKernel(const typename Dimension::Vector& x,
               const typename Dimension::SymTensor& H,
               const RKCoefficients<Dimension>& corrections) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateKernel)(*mWptr, x, H, corrections);
}

//------------------------------------------------------------------------------
// evaluateGradient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Vector
ReproducingKernel<Dimension>::
evaluateGradient(const typename Dimension::Vector& x,
                 const typename Dimension::SymTensor& H,
                 const RKCoefficients<Dimension>& corrections) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateGradient)(*mWptr, x, H, corrections);
}

//------------------------------------------------------------------------------
// evaluateHessian
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::SymTensor
ReproducingKernel<Dimension>::
evaluateHessian(const typename Dimension::Vector& x,
                const typename Dimension::SymTensor& H,
                const RKCoefficients<Dimension>& corrections) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateHessian)(*mWptr, x, H, corrections);
}

//------------------------------------------------------------------------------
// evaluateKernelAndGradient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::pair<typename Dimension::Scalar, typename Dimension::Vector>
ReproducingKernel<Dimension>::
evaluateKernelAndGradient(const typename Dimension::Vector& x,
                          const typename Dimension::SymTensor& H,
                          const RKCoefficients<Dimension>& corrections) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateKernelAndGradient)(*mWptr, x, H, corrections);
}

//------------------------------------------------------------------------------
// evaluateKernelAndGradients
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::tuple<typename Dimension::Scalar, typename Dimension::Vector, typename Dimension::Scalar>
ReproducingKernel<Dimension>::
evaluateKernelAndGradients(const typename Dimension::Vector& x,
                           const typename Dimension::SymTensor& H,
                           const RKCoefficients<Dimension>& corrections) const {
  return (*ReproducingKernelMethods<Dimension>::mEvaluateKernelAndGradients)(*mWptr, x, H, corrections);
}

//------------------------------------------------------------------------------
// computeCorrections
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
ReproducingKernel<Dimension>::
computeCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                   const FieldList<Dimension, Scalar>& volume,
                   const FieldList<Dimension, Vector>& position,
                   const FieldList<Dimension, SymTensor>& H,
                   const bool needHessian,
                   FieldList<Dimension, RKCoefficients<Dimension>>& zerothCorrections,
                   FieldList<Dimension, RKCoefficients<Dimension>>& corrections) const {
  (*ReproducingKernelMethods<Dimension>::mComputeCorrections)(connectivityMap,
                                                              *mWptr,
                                                              volume,
                                                              position,
                                                              H,
                                                              needHessian,
                                                              zerothCorrections,
                                                              corrections);
}

//------------------------------------------------------------------------------
// computeNormal
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
ReproducingKernel<Dimension>::
computeNormal(const ConnectivityMap<Dimension>& connectivityMap,
              const FieldList<Dimension, Scalar>& volume,
              const FieldList<Dimension, Vector>& position,
              const FieldList<Dimension, SymTensor>& H,
              const FieldList<Dimension, RKCoefficients<Dimension>>& corrections,
              FieldList<Dimension, Scalar>& surfaceArea,
              FieldList<Dimension, Vector>& normal) const {
  (*ReproducingKernelMethods<Dimension>::mComputeNormal)(connectivityMap,
                                                         *mWptr,
                                                         volume,
                                                         position,
                                                         H,
                                                         corrections,
                                                         surfaceArea,
                                                         normal);
}

//------------------------------------------------------------------------------
// kernel
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
ReproducingKernel<Dimension>::
kernel() const {
  return *mWptr;
}

}
