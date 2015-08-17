text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "splatMultipleFieldsMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //============================== splatFieldsMash() ==============================
    template 
    FieldListSet< Dim< %(ndim)s > >
    splatMultipleFieldsMash< Dim< %(ndim)s > >(const FieldListSet< Dim< %(ndim)s > >& fieldListSet,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& weight,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                      const TableKernel< Dim< %(ndim)s > >& kernel,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& samplePositions,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& sampleWeight,
                                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& sampleHfield,
                                      const vector<Boundary<Dim< %(ndim)s > >*>& boundaryConditions);

  }
}
"""
