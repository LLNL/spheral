namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// Access the reflection operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Tensor&
RigidBoundary<Dimension>::reflectOperator() const {
  return mReflectOperator;
}

}
}
