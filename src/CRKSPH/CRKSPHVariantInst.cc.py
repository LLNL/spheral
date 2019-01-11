text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "CRKSPHVariant.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class CRKSPHVariant< Dim< %(ndim)s > >;
}
"""
