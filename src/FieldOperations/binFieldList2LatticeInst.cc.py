text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FieldOperations/binFieldList2Lattice.cc"

namespace Spheral {

template std::vector<Dim< %(ndim)s >::Scalar> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);
template std::vector<Dim< %(ndim)s >::Vector> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);
template std::vector<Dim< %(ndim)s >::Tensor> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);
template std::vector<Dim< %(ndim)s >::SymTensor> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);

template std::vector<Dim< %(ndim)s >::Scalar> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList, const TableKernel<Dim< %(ndim)s > >& W, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);
template std::vector<Dim< %(ndim)s >::Vector> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList, const TableKernel<Dim< %(ndim)s > >& W, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);
template std::vector<Dim< %(ndim)s >::Tensor> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor>& fieldList, const TableKernel<Dim< %(ndim)s > >& W, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);
template std::vector<Dim< %(ndim)s >::SymTensor> binFieldList2Lattice<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& fieldList, const TableKernel<Dim< %(ndim)s > >& W, const Dim< %(ndim)s >::Vector& xmin, const Dim< %(ndim)s >::Vector& xmax, const std::vector<unsigned>& nsample);

}
"""
