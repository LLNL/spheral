#include "Utilities/DBC.hh"
#include <algorithm>

namespace Spheral {

VVI_IMPL_BEGIN

//------------------------------------------------------------------------------
// Construct with the given name.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::
FieldBase(typename FieldBase<Dimension>::FieldName name):
  mName(name),
  mNodeListPtr(0)
{}

//------------------------------------------------------------------------------
// Construct for a given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::
FieldBase(typename FieldBase<Dimension>::FieldName name,
          const NodeList<Dimension>& nodeList):
  mName(name),
  mNodeListPtr(&nodeList) {
#ifndef VVI_ENABLED
  mNodeListPtr->registerField(*this);
#endif
}

//------------------------------------------------------------------------------
// Copy Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::FieldBase(const FieldBase& fieldBase):
  mName(fieldBase.name()),
  mNodeListPtr(fieldBase.nodeListPtr()) {
#ifndef VVI_ENABLED
  mNodeListPtr->registerField(*this);
#endif
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
FieldBase<Dimension>::~FieldBase() {
#ifndef VVI_ENABLED
  if (mNodeListPtr) mNodeListPtr->unregisterField(*this);
#endif
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
//template<typename Dimension>
//inline
//FieldBase<Dimension>&
//FieldBase<Dimension>::operator=(const FieldBase<Dimension>& rhs) {
//  if (this != &rhs) {
//    mNodeListPtr = rhs.mNodeListPtr;
//  }
//  return *this;
//}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
FieldBase<Dimension>::operator!=(const FieldBase<Dimension>& rhs) const {
  return not (*this == rhs);
}

//------------------------------------------------------------------------------
// Get the name string.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename FieldBase<Dimension>::FieldName
FieldBase<Dimension>::name() const {
  return mName;
}

//------------------------------------------------------------------------------
// Set the name string.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::name(typename FieldBase<Dimension>::FieldName name) {
  mName = name;
}

//------------------------------------------------------------------------------
// Return a reference to this field's NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>&
FieldBase<Dimension>::nodeList() const {
  CHECK(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// Return a pointer to this field's NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>*
FieldBase<Dimension>::nodeListPtr() const {
  return mNodeListPtr;
}

//------------------------------------------------------------------------------
// Set the node list pointer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::unregisterNodeList() {
#ifndef VVI_ENABLED
  mNodeListPtr->unregisterField(*this);
#endif
  mNodeListPtr = 0;
}

//------------------------------------------------------------------------------
// Set the node list pointer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
FieldBase<Dimension>::setFieldBaseNodeList(const NodeList<Dimension>& nodeList) {
#ifndef VVI_ENABLED
  if (mNodeListPtr != 0) unregisterNodeList();
#endif
  mNodeListPtr = &nodeList;
#ifndef VVI_ENABLED
  nodeList.registerField(*this);
#endif
}


VVI_IMPL_END

} // namespace Spheral
