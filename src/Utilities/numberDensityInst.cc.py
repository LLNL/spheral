text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/numberDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

  template FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar> numberDensity<Dim< %(ndim)s > >(const DataBase<Dim< %(ndim)s > >& dataBase, const TableKernel<Dim< %(ndim)s > >& W);

}
"""
