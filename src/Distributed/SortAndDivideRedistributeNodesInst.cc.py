text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SortAndDivideRedistributeNodes.cc"

template class Spheral::PartitionSpace::SortAndDivideRedistributeNodes< Spheral::Dim< %(ndim)s > >;
"""
