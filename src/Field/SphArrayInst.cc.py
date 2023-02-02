text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Field/SphArray.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
"""

for DT in ("Dim< %(ndim)s >::Scalar",
           "Dim< %(ndim)s >::Vector",
           "Dim< %(ndim)s >::Tensor",
           "Dim< %(ndim)s >::SymTensor"):
    text += """

template class SphArrayIterator< SphArray< %(DT)s > >;
template class SphArrayIterator< SphArrayView< %(DT)s > >;
""" % {"DT" : DT}

text += """
}
"""

