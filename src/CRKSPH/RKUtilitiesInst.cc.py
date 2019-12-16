text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "RKUtilities.cc"
#include "Geometry/Dimension.hh"
"""

for order in ["ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder"]:
    text += """
namespace Spheral {
template class RKUtilities<Dim<%(ndim)s>, """
    text += """CRKOrder::%(order)s>;
}
""" % {"order" : order}
