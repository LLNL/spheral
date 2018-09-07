text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FixedSmoothingScale.cc"
#include "Geometry/Dimension.hh"

template class Spheral::FixedSmoothingScale<Spheral::Dim< %(ndim)s > >;
"""
