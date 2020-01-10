namespace Spheral {

//------------------------------------------------------------------------------
// transformationMatrix
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename ReproducingKernelMethods<Dimension>::TransformationMatrix
ReproducingKernelMethods<Dimension>::
transformationMatrix(const typename Dimension::Tensor& T,
                     const bool needHessian) const {
  TransformationMatrix result;
  (*mGetTransformationMatrix)(T, needHessian, result);
  return result;
}

//------------------------------------------------------------------------------
// applyTransformation
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
ReproducingKernelMethods<Dimension>::
applyTransformation(const typename ReproducingKernelMethods<Dimension>::TransformationMatrix& T,
                    RKCoefficients<Dimension>& corrections) const {
  (*mApplyTransformation)(T, corrections);
}

//------------------------------------------------------------------------------
// order
//------------------------------------------------------------------------------
template<typename Dimension>
inline
RKOrder
ReproducingKernelMethods<Dimension>::
order() const {
  return mOrder;
}

//------------------------------------------------------------------------------
// gradCorrectionsSize
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
ReproducingKernelMethods<Dimension>::
gradCorrectionsSize() const {
  return mGradCorrectionsSize;
}

//------------------------------------------------------------------------------
// hessCorrectionsSize
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
ReproducingKernelMethods<Dimension>::
hessCorrectionsSize() const {
  return mHessCorrectionsSize;
}

}
