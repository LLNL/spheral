text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SortAndDivideRedistributeNodes.cc"

template class Spheral::SortAndDivideRedistributeNodes< Spheral::Dim< %(ndim)s > >;
"""
