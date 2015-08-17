text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "divergencePairWise.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //============================== gradient() ==============================
    template 
    FieldList<Dim< %(ndim)s >, MathTraits<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>::DivergenceType> 
    divergencePairWise<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& density,
                                               const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                               const TableKernel< Dim< %(ndim)s > >& kernel);

  }
}
"""
