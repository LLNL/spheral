text = """
//------------------------------------------------------------------------------
// Initialize the static instance pointer.
//------------------------------------------------------------------------------
#include "NodeListRegistrar.cc"

template<typename Dimension> Spheral::NodeListRegistrar<Dimension>* Spheral::NodeListRegistrar<Dimension>::mInstancePtr = 0;
template class Spheral::NodeListRegistrar<Spheral::Dim< %(ndim)s > >;
"""
