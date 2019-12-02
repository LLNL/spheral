text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SuperiorRKUtilities.cc"
#include "Geometry/Dimension.hh"
"""

for order in ["ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder"]:
    text += """
namespace Spheral {
template class SuperiorRKUtilities<Dim<%(ndim)s>, """
    text += """CRKOrder::%(order)s>;
}
""" % {"order" : order}
