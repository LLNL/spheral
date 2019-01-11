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

}
