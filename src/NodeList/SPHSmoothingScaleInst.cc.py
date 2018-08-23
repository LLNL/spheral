text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/SPHSmoothingScale.cc"
#include "Geometry/Dimension.hh"

template class Spheral::NodeSpace::SPHSmoothingScale<Spheral::Dim< %(ndim)s > >;
"""
