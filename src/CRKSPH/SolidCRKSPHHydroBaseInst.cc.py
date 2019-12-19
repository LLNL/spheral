text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SolidCRKSPHHydroBase.cc"

namespace Spheral {
template class SolidCRKSPHHydroBase< Dim< %(ndim)s > >;
}
"""
