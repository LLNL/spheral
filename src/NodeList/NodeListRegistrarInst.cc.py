text = """
//------------------------------------------------------------------------------
// Initialize the static instance pointer.
//------------------------------------------------------------------------------
#include "NodeList/NodeListRegistrar.cc"
template class Spheral::NodeListRegistrar<Spheral::Dim<%(ndim)s>>;
"""
