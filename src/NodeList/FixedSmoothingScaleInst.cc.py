text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/FixedSmoothingScale.cc"
#include "Geometry/Dimension.hh"

template class Spheral::NodeSpace::FixedSmoothingScale<Spheral::Dim< %(ndim)s > >;
"""
