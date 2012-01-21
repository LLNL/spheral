#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Neighbor/GridCellIndex.hh"

// BPL includes
#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;

//------------------------------------------------------------------------------
// Wrap std::vector of primitive types.
//------------------------------------------------------------------------------
void
wrapSTLVectorPrimitives() {

  // Single level vectors.
  class_<std::vector<char> >("vector_of_char")
    .def(boost::python::vector_indexing_suite<std::vector<char> >())
    ;

  class_<std::vector<int> >("vector_of_int")
    .def(boost::python::vector_indexing_suite<std::vector<int> >())
    ;

  class_<std::vector<bool> >("vector_of_bool")
    .def(boost::python::vector_indexing_suite<std::vector<bool> >())
    ;

  class_<std::vector<float> >("vector_of_float")
    .def(boost::python::vector_indexing_suite<std::vector<float> >())
    ;

  class_<std::vector<double> >("vector_of_double")
    .def(boost::python::vector_indexing_suite<std::vector<double> >())
    ;

  class_<std::vector<std::string> >("vector_of_string")
    .def(boost::python::vector_indexing_suite<std::vector<std::string> >())
    ;

  class_<std::vector<Dim<1>::Vector> >("vector_of_Vector1d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<1>::Vector> >())
    ;

  class_<std::vector<Dim<1>::Tensor> >("vector_of_Tensor1d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<1>::Tensor> >())
    ;

  class_<std::vector<Dim<1>::SymTensor> >("vector_of_SymTensor1d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<1>::SymTensor> >())
    ;

  class_<std::vector<Dim<1>::ThirdRankTensor> >("vector_of_ThirdRankTensor1d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<1>::ThirdRankTensor> >())
    ;

  class_<std::vector<Dim<2>::Vector> >("vector_of_Vector2d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<2>::Vector> >())
    ;

  class_<std::vector<Dim<2>::Tensor> >("vector_of_Tensor2d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<2>::Tensor> >())
    ;

  class_<std::vector<Dim<2>::SymTensor> >("vector_of_SymTensor2d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<2>::SymTensor> >())
    ;

  class_<std::vector<Dim<2>::ThirdRankTensor> >("vector_of_ThirdRankTensor2d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<2>::ThirdRankTensor> >())
    ;

  class_<std::vector<Dim<3>::Vector> >("vector_of_Vector3d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<3>::Vector> >())
    ;

  class_<std::vector<Dim<3>::Tensor> >("vector_of_Tensor3d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<3>::Tensor> >())
    ;

  class_<std::vector<Dim<3>::SymTensor> >("vector_of_SymTensor3d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<3>::SymTensor> >())
    ;

  class_<std::vector<Dim<3>::ThirdRankTensor> >("vector_of_ThirdRankTensor3d")
    .def(boost::python::vector_indexing_suite<std::vector<Dim<3>::ThirdRankTensor> >())
    ;

  class_<std::vector<GeomPlane<Dim<1> > > >("vector_of_Plane1d")
    .def(boost::python::vector_indexing_suite<std::vector<GeomPlane<Dim<1> > > >())
    ;

  class_<std::vector<GeomPlane<Dim<2> > > >("vector_of_Plane2d")
    .def(boost::python::vector_indexing_suite<std::vector<GeomPlane<Dim<2> > > >())
    ;

  class_<std::vector<GeomPlane<Dim<3> > > >("vector_of_Plane3d")
    .def(boost::python::vector_indexing_suite<std::vector<GeomPlane<Dim<3> > > >())
    ;

  class_<std::vector<NeighborSpace::GridCellIndex<Dim<1> > > >("vector_of_GridCellIndex1d")
    .def(boost::python::vector_indexing_suite<std::vector< NeighborSpace::GridCellIndex<Dim<1> > > >())
    ;

  class_<std::vector<NeighborSpace::GridCellIndex<Dim<2> > > >("vector_of_GridCellIndex2d")
    .def(boost::python::vector_indexing_suite<std::vector< NeighborSpace::GridCellIndex<Dim<2> > > >())
    ;

  class_<std::vector<NeighborSpace::GridCellIndex<Dim<3> > > >("vector_of_GridCellIndex3d")
    .def(boost::python::vector_indexing_suite<std::vector< NeighborSpace::GridCellIndex<Dim<3> > > >())
    ;

  // Two level vectors.
  class_<std::vector< std::vector<int> > >("vector_of_vector_of_int")
    .def(boost::python::vector_indexing_suite<std::vector< std::vector<int> > >())
    ;

  class_<std::vector<std::vector<double> > >("vector_of_vector_of_double")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<double> > >())
    ;

  class_<std::vector<std::vector<Dim<1>::Vector> > >("vector_of_vector_of_Vector1d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<1>::Vector> > >())
    ;

  class_<std::vector<std::vector<Dim<1>::Tensor> > >("vector_of_vector_of_Tensor1d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<1>::Tensor> > >())
    ;

  class_<std::vector<std::vector<Dim<1>::SymTensor> > >("vector_of_vector_of_SymTensor1d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<1>::SymTensor> > >())
    ;

  class_<std::vector<std::vector<Dim<1>::ThirdRankTensor> > >("vector_of_vector_of_ThirdRankTensor1d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<1>::ThirdRankTensor> > >())
    ;

  class_<std::vector<std::vector<Dim<2>::Vector> > >("vector_of_vector_of_Vector2d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<2>::Vector> > >())
    ;

  class_<std::vector<std::vector<Dim<2>::Tensor> > >("vector_of_vector_of_Tensor2d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<2>::Tensor> > >())
    ;

  class_<std::vector<std::vector<Dim<2>::SymTensor> > >("vector_of_vector_of_SymTensor2d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<2>::SymTensor> > >())
    ;

  class_<std::vector<std::vector<Dim<2>::ThirdRankTensor> > >("vector_of_vector_of_ThirdRankTensor2d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<2>::ThirdRankTensor> > >())
    ;

  class_<std::vector<std::vector<Dim<3>::Vector> > >("vector_of_vector_of_Vector3d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<3>::Vector> > >())
    ;

  class_<std::vector<std::vector<Dim<3>::Tensor> > >("vector_of_vector_of_Tensor3d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<3>::Tensor> > >())
    ;

  class_<std::vector<std::vector<Dim<3>::SymTensor> > >("vector_of_vector_of_SymTensor3d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<3>::SymTensor> > >())
    ;

  class_<std::vector<std::vector<Dim<3>::ThirdRankTensor> > >("vector_of_vector_of_ThirdRankTensor3d")
    .def(boost::python::vector_indexing_suite<std::vector<std::vector<Dim<3>::ThirdRankTensor> > >())
    ;

  // Three level vectors.
  class_<std::vector< std::vector< std::vector<int> > > >("vector_of_vector_of_vector_of_int")
    .def(boost::python::vector_indexing_suite<std::vector< std::vector< std::vector<int> > > >())
    ;

}


}
