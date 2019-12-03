text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "SuperiorRKCorrections.cc"
#include "Geometry/Dimension.hh"
"""

for order in ["ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder"]:
    text += """
namespace Spheral {
template class SuperiorRKCorrections<Dim<%(ndim)s>, """
    text += """CRKOrder::%(order)s>;
}
""" % {"order" : order}
