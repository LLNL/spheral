text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "RK/RKCorrections.cc"
#include "Geometry/Dimension.hh"
"""

for order in ["ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder"]:
    text += """
namespace Spheral {
template class RKCorrections<Dim<%(ndim)s>, """
    text += """RKOrder::%(order)s>;
}
""" % {"order" : order}
