#ifndef __PBGWRAP_REFERENCE__
#define __PBGWRAP_REFERENCE__

#include <functional>

namespace Spheral {

//------------------------------------------------------------------------------
// We use these methods to get values that Spheral wants to serve up as 
// references as pointers, since pybindgen handles pointers without copies but 
// not references.  I know this is clunky!
//------------------------------------------------------------------------------
template<typename Object, typename ReturnType, const ReturnType& (Object::* AccessMethod)() const>
inline
ReturnType*
const_reference_as_pointer(const Object& self) {
  return const_cast<ReturnType*>(&((self.*AccessMethod)()));
}

template<typename Object, typename ReturnType, ReturnType& (Object::* AccessMethod)()>
inline
ReturnType*
reference_as_pointer(Object& self) {
  return &((self.*AccessMethod)());
}

template<typename Object, typename ReturnType, ReturnType& (Object::* AccessMethod)() const>
inline
ReturnType*
reference_as_pointer(Object& self) {
  return &((self.*AccessMethod)());
}

}

#endif
