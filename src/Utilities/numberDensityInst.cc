//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Utilities/numberDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

  template FieldList<Dim<1>, Dim<1>::Scalar> numberDensity<Dim<1> >(const DataBase<Dim<1> >& dataBase, const TableKernel<Dim<1> >& W);

#endif

#if defined(SPHERAL_ENABLE_2D)

  template FieldList<Dim<2>, Dim<2>::Scalar> numberDensity<Dim<2> >(const DataBase<Dim<2> >& dataBase, const TableKernel<Dim<2> >& W);

#endif

#if defined(SPHERAL_ENABLE_3D)

  template FieldList<Dim<3>, Dim<3>::Scalar> numberDensity<Dim<3> >(const DataBase<Dim<3> >& dataBase, const TableKernel<Dim<3> >& W);

#endif
}