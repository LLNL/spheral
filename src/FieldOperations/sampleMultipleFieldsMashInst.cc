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

    template 
    FieldListSet< Dim<2> >
    sampleMultipleFieldsMash< Dim<2> >(const FieldListSet< Dim<2> >& fieldListSet,
                                       const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                       const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                       const TableKernel< Dim<2> >& kernel,
                                       const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                       const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                       const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield);

    template 
    FieldListSet< Dim<3> >
    sampleMultipleFieldsMash< Dim<3> >(const FieldListSet< Dim<3> >& fieldListSet,
                                       const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                       const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                       const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                       const TableKernel< Dim<3> >& kernel,
                                       const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                       const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                       const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield);


  }
}
