//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/computeShepardsInterpolation.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
"""

for DT in ("Dim<1>::Scalar",
           "Dim<1>::Vector",
           "Dim<1>::Tensor",
           "Dim<1>::SymTensor"):
    text += """

template
FieldList<Dim< %1 >, %(DT)s>
computeShepardsInterpolation(const FieldList<Dim< %1 >, %(DT)s>& fieldList,
                             const ConnectivityMap<Dim< %1 > >&,
                             const TableKernel<Dim< %1 > >&,
                             const FieldList<Dim< %1 >, Dim< %1 >::Vector>& position,
                             const FieldList<Dim< %1 >, Dim< %1 >::SymTensor>& H,
                             const FieldList<Dim< %1 >, Dim< %1 >::Scalar>& weight);

""" % {"DT" : DT}

text += """
#endif

#if defined(SPHERAL_ENABLE_2D)
"""

for DT in ("Dim<2>::Scalar",
           "Dim<2>::Vector",
           "Dim<2>::Tensor",
           "Dim<2>::SymTensor"):
    text += """

template
FieldList<Dim< %2 >, %(DT)s>
computeShepardsInterpolation(const FieldList<Dim< %2 >, %(DT)s>& fieldList,
                             const ConnectivityMap<Dim< %2 > >&,
                             const TableKernel<Dim< %2 > >&,
                             const FieldList<Dim< %2 >, Dim< %2 >::Vector>& position,
                             const FieldList<Dim< %2 >, Dim< %2 >::SymTensor>& H,
                             const FieldList<Dim< %2 >, Dim< %2 >::Scalar>& weight);

""" % {"DT" : DT}

text += """
#endif

#if defined(SPHERAL_ENABLE_3D)
"""

for DT in ("Dim<3>::Scalar",
           "Dim<3>::Vector",
           "Dim<3>::Tensor",
           "Dim<3>::SymTensor"):
    text += """

template
FieldList<Dim< %3 >, %(DT)s>
computeShepardsInterpolation(const FieldList<Dim< %3 >, %(DT)s>& fieldList,
                             const ConnectivityMap<Dim< %3 > >&,
                             const TableKernel<Dim< %3 > >&,
                             const FieldList<Dim< %3 >, Dim< %3 >::Vector>& position,
                             const FieldList<Dim< %3 >, Dim< %3 >::SymTensor>& H,
                             const FieldList<Dim< %3 >, Dim< %3 >::Scalar>& weight);

""" % {"DT" : DT}

text += """
#endif
}