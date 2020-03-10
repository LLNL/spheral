text = """
//------------------------------------------------------------------------------
// Initialize the static instance pointer.
//------------------------------------------------------------------------------
#include "NodeList/NodeListRegistrar.cc"

template<> Spheral::NodeListRegistrar<Spheral::Dim<%(ndim)s>>* Spheral::NodeListRegistrar<Spheral::Dim<%(ndim)s>>::mInstancePtr = 0;
template class Spheral::NodeListRegistrar<Spheral::Dim<%(ndim)s>>;
"""
