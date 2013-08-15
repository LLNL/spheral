//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "numberDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

  template FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar> numberDensity<Dim<1> >(const DataBaseSpace::DataBase<Dim<1> >& dataBase, const KernelSpace::TableKernel<Dim<1> >& W);
  template FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar> numberDensity<Dim<2> >(const DataBaseSpace::DataBase<Dim<2> >& dataBase, const KernelSpace::TableKernel<Dim<2> >& W);
  template FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar> numberDensity<Dim<3> >(const DataBaseSpace::DataBase<Dim<3> >& dataBase, const KernelSpace::TableKernel<Dim<3> >& W);

}
