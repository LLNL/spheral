//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "FieldOperations/gradDivVectorFieldListGolden2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

//========================== gradDivVectorFieldListGolden2() ==========================
template
FieldList<Dim<1>, Dim<1>::Vector>
gradDivVectorFieldListGolden2< Dim<1> >
(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
 const FieldList<Dim<1>, Dim<1>::Vector>& position,
 const FieldList<Dim<1>, Dim<1>::Scalar>& weight,
 const FieldList<Dim<1>, Dim<1>::Scalar>& mass,
 const FieldList<Dim<1>, Dim<1>::Scalar>& rho,
 const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
 const TableKernel< Dim<1> >& kernel);

#endif

#if defined(SPHERAL_ENABLE_2D)

//========================== gradDivVectorFieldListGolden2() ==========================
template
FieldList<Dim<2>, Dim<2>::Vector>
gradDivVectorFieldListGolden2< Dim<2> >
(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
 const FieldList<Dim<2>, Dim<2>::Vector>& position,
 const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
 const FieldList<Dim<2>, Dim<2>::Scalar>& mass,
 const FieldList<Dim<2>, Dim<2>::Scalar>& rho,
 const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
 const TableKernel< Dim<2> >& kernel);

#endif

#if defined(SPHERAL_ENABLE_3D)

//========================== gradDivVectorFieldListGolden2() ==========================
template
FieldList<Dim<3>, Dim<3>::Vector>
gradDivVectorFieldListGolden2< Dim<3> >
(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
 const FieldList<Dim<3>, Dim<3>::Vector>& position,
 const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
 const FieldList<Dim<3>, Dim<3>::Scalar>& mass,
 const FieldList<Dim<3>, Dim<3>::Scalar>& rho,
 const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
 const TableKernel< Dim<3> >& kernel);

#endif
}