text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SPHSmoothingScale.cc"
#include "Geometry/Dimension.hh"

template class Spheral::SPHSmoothingScale<Spheral::Dim< %(ndim)s > >;
"""
