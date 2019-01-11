text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/SPHSmoothingScale.cc"
#include "Geometry/Dimension.hh"

template class Spheral::SPHSmoothingScale<Spheral::Dim< %(ndim)s > >;
"""
