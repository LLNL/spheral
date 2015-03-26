//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPHSmoothingScale.cc"
#include "Geometry/Dimension.hh"

template class Spheral::NodeSpace::SPHSmoothingScale<Spheral::Dim<1> >;
template class Spheral::NodeSpace::SPHSmoothingScale<Spheral::Dim<2> >;
template class Spheral::NodeSpace::SPHSmoothingScale<Spheral::Dim<3> >;
