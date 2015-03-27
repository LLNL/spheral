//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "binFieldList2Lattice.cc"

namespace Spheral {
  namespace FieldSpace {

    template std::vector<Dim<1>::Scalar> binFieldList2Lattice<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<1>::Vector> binFieldList2Lattice<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<1>::Tensor> binFieldList2Lattice<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<1>::SymTensor> binFieldList2Lattice<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);

    template std::vector<Dim<2>::Scalar> binFieldList2Lattice<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<2>::Vector> binFieldList2Lattice<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<2>::Tensor> binFieldList2Lattice<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<2>::SymTensor> binFieldList2Lattice<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);

    template std::vector<Dim<3>::Scalar> binFieldList2Lattice<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<3>::Vector> binFieldList2Lattice<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<3>::Tensor> binFieldList2Lattice<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<3>::SymTensor> binFieldList2Lattice<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);

    template std::vector<Dim<1>::Scalar> binFieldList2Lattice<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList, const TableKernel<Dim<1> >& W, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<1>::Vector> binFieldList2Lattice<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList, const TableKernel<Dim<1> >& W, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<1>::Tensor> binFieldList2Lattice<Dim<1>, Dim<1>::Tensor>(const FieldList<Dim<1>, Dim<1>::Tensor>& fieldList, const TableKernel<Dim<1> >& W, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<1>::SymTensor> binFieldList2Lattice<Dim<1>, Dim<1>::SymTensor>(const FieldList<Dim<1>, Dim<1>::SymTensor>& fieldList, const TableKernel<Dim<1> >& W, const Dim<1>::Vector& xmin, const Dim<1>::Vector& xmax, const std::vector<unsigned>& nsample);

    template std::vector<Dim<2>::Scalar> binFieldList2Lattice<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList, const TableKernel<Dim<2> >& W, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<2>::Vector> binFieldList2Lattice<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList, const TableKernel<Dim<2> >& W, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<2>::Tensor> binFieldList2Lattice<Dim<2>, Dim<2>::Tensor>(const FieldList<Dim<2>, Dim<2>::Tensor>& fieldList, const TableKernel<Dim<2> >& W, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<2>::SymTensor> binFieldList2Lattice<Dim<2>, Dim<2>::SymTensor>(const FieldList<Dim<2>, Dim<2>::SymTensor>& fieldList, const TableKernel<Dim<2> >& W, const Dim<2>::Vector& xmin, const Dim<2>::Vector& xmax, const std::vector<unsigned>& nsample);

    template std::vector<Dim<3>::Scalar> binFieldList2Lattice<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList, const TableKernel<Dim<3> >& W, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<3>::Vector> binFieldList2Lattice<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList, const TableKernel<Dim<3> >& W, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<3>::Tensor> binFieldList2Lattice<Dim<3>, Dim<3>::Tensor>(const FieldList<Dim<3>, Dim<3>::Tensor>& fieldList, const TableKernel<Dim<3> >& W, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);
    template std::vector<Dim<3>::SymTensor> binFieldList2Lattice<Dim<3>, Dim<3>::SymTensor>(const FieldList<Dim<3>, Dim<3>::SymTensor>& fieldList, const TableKernel<Dim<3> >& W, const Dim<3>::Vector& xmin, const Dim<3>::Vector& xmax, const std::vector<unsigned>& nsample);

  }
}
