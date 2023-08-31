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

#for FV in ("FieldView< Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar >",
#           "FieldView< Dim< %(ndim)s >, Dim< %(ndim)s >::Vector >",
#           "FieldView< Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor >",
#           "FieldView< Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor >"):
#    text += """
#template class SphArrayFieldIterator< SphArray< %(FV)s > >;
#template class SphArrayFieldIterator< SphArrayView< %(FV)s > >;
#""" % {"FV" : FV}

text += """
}
"""

