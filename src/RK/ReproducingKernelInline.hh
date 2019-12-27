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
  return (*mEvaluateBaseKernel)(mW, x, H);
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
  return (*mEvaluateBaseGradient)(mW, x, H);
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
  return (*mEvaluateBaseHessian)(mW, x, H);
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
  return (*mEvaluateBaseKernelAndGradient)(mW, x, H);
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
               const std::vector<double>& corrections) const {
  return (*mEvaluateKernel)(mW, x, H, corrections);
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
                 const std::vector<double>& corrections) const {
  return (*mEvaluateGradient)(mW, x, H, corrections);
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
                const std::vector<double>& corrections) const {
  return (*mEvaluateHessian)(mW, x, H, corrections);
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
                          const std::vector<double>& corrections) const {
  return (*mEvaluateKernelAndGradient)(mW, x, H, corrections);
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
                           const std::vector<double>& corrections) const {
  return (*mEvaluateKernelAndGradients)(mW, x, H, corrections);
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
                   FieldList<Dimension, std::vector<double>>& zerothCorrections,
                   FieldList<Dimension, std::vector<double>>& corrections) {
  (*mComputeCorrections)(connectivityMap,
                         mW,
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
              const FieldList<Dimension, std::vector<double>>& corrections,
              FieldList<Dimension, Scalar>& surfaceArea,
              FieldList<Dimension, Vector>& normal) {
  (*mComputeNormal)(connectivityMap,
                    mW,
                    volume,
                    position,
                    H,
                    corrections,
                    surfaceArea,
                    normal);
}

}
