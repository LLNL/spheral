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
    FieldListSet< Dim<1> >
    splatMultipleFieldsMash< Dim<1> >(const FieldListSet< Dim<1> >& fieldListSet,
                                      const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                      const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                                      const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                      const TableKernel< Dim<1> >& kernel,
                                      const FieldList<Dim<1>, Dim<1>::Vector>& samplePositions,
                                      const FieldList<Dim<1>, Dim<1>::Scalar>& sampleWeight,
                                      const FieldList<Dim<1>, Dim<1>::SymTensor>& sampleHfield,
                                      const vector<Boundary<Dim<1> >*>& boundaryConditions);

    template 
    FieldListSet< Dim<2> >
    splatMultipleFieldsMash< Dim<2> >(const FieldListSet< Dim<2> >& fieldListSet,
                                      const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                      const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                                      const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                      const TableKernel< Dim<2> >& kernel,
                                      const FieldList<Dim<2>, Dim<2>::Vector>& samplePositions,
                                      const FieldList<Dim<2>, Dim<2>::Scalar>& sampleWeight,
                                      const FieldList<Dim<2>, Dim<2>::SymTensor>& sampleHfield,
                                      const vector<Boundary<Dim<2> >*>& boundaryConditions);

    template 
    FieldListSet< Dim<3> >
    splatMultipleFieldsMash< Dim<3> >(const FieldListSet< Dim<3> >& fieldListSet,
                                      const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                      const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                                      const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                      const TableKernel< Dim<3> >& kernel,
                                      const FieldList<Dim<3>, Dim<3>::Vector>& samplePositions,
                                      const FieldList<Dim<3>, Dim<3>::Scalar>& sampleWeight,
                                      const FieldList<Dim<3>, Dim<3>::SymTensor>& sampleHfield,
                                      const vector<Boundary<Dim<3> >*>& boundaryConditions);

  }
}
