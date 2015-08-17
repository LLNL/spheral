text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "gradDivVectorFieldListGolden.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //========================== gradDivVectorFieldListGolden() ==========================
    template 
    FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
    gradDivVectorFieldListGolden< Dim< %(ndim)s > >
    (const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& mass,
     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& rho,
     const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
     const TableKernel< Dim< %(ndim)s > >& kernel);

  }
}
"""
