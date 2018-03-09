text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "numberDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

  template FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> numberDensity<Dim< %(ndim)s > >(const DataBaseSpace::DataBase<Dim< %(ndim)s > >& dataBase, const KernelSpace::TableKernel<Dim< %(ndim)s > >& W);

}
"""
