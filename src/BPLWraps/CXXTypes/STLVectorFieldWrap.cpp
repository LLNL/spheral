#include "Field/Field.hh"
#include "Geometry/Dimension.hh"

// BPL includes
#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;
using namespace Spheral::FieldSpace;

//------------------------------------------------------------------------------
// Wrap std::vector of Field.
//------------------------------------------------------------------------------
void
wrapSTLVectorField() {

  //----------------------------------------------------------------------------
  // 1-D
  class_<std::vector<Field<Dim<1>, Dim<1>::Scalar> > >("vector_of_ScalarField1d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<1>, Dim<1>::Scalar> > >())
    ;

  class_<std::vector<Field<Dim<1>, Dim<1>::Vector> > >("vector_of_VectorField1d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<1>, Dim<1>::Vector> > >())
    ;

  class_<std::vector<Field<Dim<1>, Dim<1>::Tensor> > >("vector_of_TensorField1d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<1>, Dim<1>::Tensor> > >())
    ;

  class_<std::vector<Field<Dim<1>, Dim<1>::SymTensor> > >("vector_of_SymTensorField1d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<1>, Dim<1>::SymTensor> > >())
    ;

  //----------------------------------------------------------------------------
  // 2-D
  class_<std::vector<Field<Dim<2>, Dim<2>::Scalar> > >("vector_of_ScalarField2d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<2>, Dim<2>::Scalar> > >())
    ;

  class_<std::vector<Field<Dim<2>, Dim<2>::Vector> > >("vector_of_VectorField2d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<2>, Dim<2>::Vector> > >())
    ;

  class_<std::vector<Field<Dim<2>, Dim<2>::Tensor> > >("vector_of_TensorField2d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<2>, Dim<2>::Tensor> > >())
    ;

  class_<std::vector<Field<Dim<2>, Dim<2>::SymTensor> > >("vector_of_SymTensorField2d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<2>, Dim<2>::SymTensor> > >())
    ;

  //----------------------------------------------------------------------------
  // 3-D
  class_<std::vector<Field<Dim<3>, Dim<3>::Scalar> > >("vector_of_ScalarField3d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<3>, Dim<3>::Scalar> > >())
    ;

  class_<std::vector<Field<Dim<3>, Dim<3>::Vector> > >("vector_of_VectorField3d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<3>, Dim<3>::Vector> > >())
    ;

  class_<std::vector<Field<Dim<3>, Dim<3>::Tensor> > >("vector_of_TensorField3d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<3>, Dim<3>::Tensor> > >())
    ;

  class_<std::vector<Field<Dim<3>, Dim<3>::SymTensor> > >("vector_of_SymTensorField3d")
    .def(boost::python::vector_indexing_suite<std::vector<Field<Dim<3>, Dim<3>::SymTensor> > >())
    ;

}

}
