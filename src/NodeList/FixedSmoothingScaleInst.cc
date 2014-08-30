//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FixedSmoothingScale.cc"
#include "Geometry/Dimension.hh"

template class Spheral::NodeSpace::FixedSmoothingScale<Spheral::Dim<1> >;
template class Spheral::NodeSpace::FixedSmoothingScale<Spheral::Dim<2> >;
template class Spheral::NodeSpace::FixedSmoothingScale<Spheral::Dim<3> >;
