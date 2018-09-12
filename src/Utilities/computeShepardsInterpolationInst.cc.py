text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/computeShepardsInterpolation.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
"""

for DT in ("Dim< %(ndim)s >::Scalar",
           "Dim< %(ndim)s >::Vector",
           "Dim< %(ndim)s >::Tensor",
           "Dim< %(ndim)s >::SymTensor"):
    text += """

template
FieldList<Dim< %%(ndim)s >, %(DT)s>
computeShepardsInterpolation(const FieldList<Dim< %%(ndim)s >, %(DT)s>& fieldList,
                             const ConnectivityMap<Dim< %%(ndim)s > >&, 
                             const TableKernel<Dim< %%(ndim)s > >&, 
                             const FieldList<Dim< %%(ndim)s >, Dim< %%(ndim)s >::Vector>& position,
                             const FieldList<Dim< %%(ndim)s >, Dim< %%(ndim)s >::SymTensor>& H,
                             const FieldList<Dim< %%(ndim)s >, Dim< %%(ndim)s >::Scalar>& weight);

""" % {"DT" : DT}

text += """
}
"""
