text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "RK/RKUtilities.cc"
#include "Geometry/Dimension.hh"
"""

for order in ["CubicOrder"]:
    text += """
namespace Spheral {
template class RKUtilities<Dim<%(ndim)s>, """
    text += """RKOrder::%(order)s>;
}
""" % {"order" : order}
