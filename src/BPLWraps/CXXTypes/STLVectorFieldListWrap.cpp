#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"

// BPL includes
#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;
using namespace Spheral::FieldSpace;

//------------------------------------------------------------------------------
// Wrap std::vector of FieldList.
//------------------------------------------------------------------------------
void
wrapSTLVectorFieldList() {

  //----------------------------------------------------------------------------
  // 1-D
  class_<std::vector<FieldList<Dim<1>, Dim<1>::Scalar> > >("vector_of_ScalarFieldList1d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<1>, Dim<1>::Scalar> > >())
    ;

  class_<std::vector<FieldList<Dim<1>, Dim<1>::Vector> > >("vector_of_VectorFieldList1d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<1>, Dim<1>::Vector> > >())
    ;

  class_<std::vector<FieldList<Dim<1>, Dim<1>::Tensor> > >("vector_of_TensorFieldList1d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<1>, Dim<1>::Tensor> > >())
    ;

  class_<std::vector<FieldList<Dim<1>, Dim<1>::SymTensor> > >("vector_of_SymTensorFieldList1d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<1>, Dim<1>::SymTensor> > >())
    ;

  //----------------------------------------------------------------------------
  // 2-D
  class_<std::vector<FieldList<Dim<2>, Dim<2>::Scalar> > >("vector_of_ScalarFieldList2d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<2>, Dim<2>::Scalar> > >())
    ;

  class_<std::vector<FieldList<Dim<2>, Dim<2>::Vector> > >("vector_of_VectorFieldList2d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<2>, Dim<2>::Vector> > >())
    ;

  class_<std::vector<FieldList<Dim<2>, Dim<2>::Tensor> > >("vector_of_TensorFieldList2d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<2>, Dim<2>::Tensor> > >())
    ;

  class_<std::vector<FieldList<Dim<2>, Dim<2>::SymTensor> > >("vector_of_SymTensorFieldList2d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<2>, Dim<2>::SymTensor> > >())
    ;

  //----------------------------------------------------------------------------
  // 3-D
  class_<std::vector<FieldList<Dim<3>, Dim<3>::Scalar> > >("vector_of_ScalarFieldList3d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<3>, Dim<3>::Scalar> > >())
    ;

  class_<std::vector<FieldList<Dim<3>, Dim<3>::Vector> > >("vector_of_VectorFieldList3d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<3>, Dim<3>::Vector> > >())
    ;

  class_<std::vector<FieldList<Dim<3>, Dim<3>::Tensor> > >("vector_of_TensorFieldList3d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<3>, Dim<3>::Tensor> > >())
    ;

  class_<std::vector<FieldList<Dim<3>, Dim<3>::SymTensor> > >("vector_of_SymTensorFieldList3d")
    .def(boost::python::vector_indexing_suite<std::vector<FieldList<Dim<3>, Dim<3>::SymTensor> > >())
    ;

}

}
