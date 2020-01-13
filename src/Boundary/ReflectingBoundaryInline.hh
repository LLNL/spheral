namespace Spheral {

//------------------------------------------------------------------------------
// Access the reflection operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Tensor&
ReflectingBoundary<Dimension>::reflectOperator() const {
  return mReflectOperator;
}

//------------------------------------------------------------------------------
// Access the RK coefficients reflection operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename ReflectingBoundary<Dimension>::TransformationMatrix&
ReflectingBoundary<Dimension>::rkReflectOperator(const RKOrder order,
                                                 const bool useHessian) const {
  const auto itr = mrkReflectOperators.find(order);
  CHECK(itr != mrkReflectOperators.end());
  return useHessian ? itr->second.second : itr->second.first;
}

}
