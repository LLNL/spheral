text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "sampleMultipleFieldsMash.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace FieldSpace {

    using KernelSpace::TableKernel;

    //============================== sampleFieldsMash() ==============================
    template 
    FieldListSet< Dim<1> >
    sampleMultipleFieldsMash< Dim<1> >(const FieldListSet< Dim<1> >& fieldListSet,
                                       const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                       const TableKernel< Dim<1> >& kernel,
                                       const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                       const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                       const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield);

  }
}
"""
