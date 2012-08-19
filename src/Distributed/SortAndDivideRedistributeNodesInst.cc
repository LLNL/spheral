//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SortAndDivideRedistributeNodes.cc"

template class Spheral::PartitionSpace::SortAndDivideRedistributeNodes< Dim<1> >;
template class Spheral::PartitionSpace::SortAndDivideRedistributeNodes< Dim<2> >;
template class Spheral::PartitionSpace::SortAndDivideRedistributeNodes< Dim<3> >;
