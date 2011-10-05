#ifndef __GCCXML__
#include "boost/python/class.hpp"
#else
#include "fakeboost.hh"
#endif

//------------------------------------------------------------------------------
// Intermediary to insulate us (somewhat) from the BPL changing interface for
// attribute wraps.
//------------------------------------------------------------------------------
namespace Spheral {

// Version where we need to explictly state a return policy.
template<typename FunctionPtr, typename ReturnPolicy>
boost::python::object
wrapAttr(FunctionPtr funcPtr, ReturnPolicy policy) {
  return boost::python::object((boost::python::make_function(funcPtr, policy)));
}

// Specialized for the case where no return policy is needed.
template<typename FunctionPtr>
boost::python::object
wrapAttr(FunctionPtr funcPtr) {
  return boost::python::object((boost::python::make_function(funcPtr)));
}

// // Version where we need to explictly state a return policy.
// template<typename FunctionPtr, typename ReturnPolicy>
// boost::python::object
// wrapAttr(FunctionPtr funcPtr, ReturnPolicy policy) {
//   return boost::python::object(boost::python::handle<>(boost::python::make_function(funcPtr, policy)));
// }

// // Specialized for the case where no return policy is needed.
// template<typename FunctionPtr>
// boost::python::object
// wrapAttr(FunctionPtr funcPtr) {
//   return boost::python::object(boost::python::handle<>(boost::python::make_function(funcPtr)));
// }

}
