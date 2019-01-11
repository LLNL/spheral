text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/SortAndDivideRedistributeNodes.cc"

template class Spheral::SortAndDivideRedistributeNodes< Spheral::Dim< %(ndim)s > >;
"""
