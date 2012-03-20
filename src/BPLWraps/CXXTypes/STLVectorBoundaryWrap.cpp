#include "Boundary/Boundary.hh"
#include "Geometry/Dimension.hh"

// BPL includes
#include "boost/python.hpp"
#include "boost/python/iterator.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;
using namespace Spheral::BoundarySpace;

//------------------------------------------------------------------------------
// Wrap std::vector of Boundary.
//------------------------------------------------------------------------------
void
wrapSTLVectorBoundary() {

//   class_<std::vector<Boundary<Dim<1> >*> >("vector_of_Boundary1d_ptr")
//     .def("__iter__", boost::python::iterator<std::vector<Boundary<Dim<1> >*>,
//          return_value_policy<reference_existing_object> >())
//     ;

  class_<std::vector<Boundary<Dim<1> >*> >("vector_of_Boundary1d_ptr")
    .def(boost::python::vector_indexing_suite<std::vector<Boundary<Dim<1> >*> >())
    .def("__iter__", boost::python::iterator<std::vector<Boundary<Dim<1> >*>,
         return_internal_reference<> >())
    ;

  class_<std::vector<Boundary<Dim<2> >*> >("vector_of_Boundary2d_ptr")
    .def(boost::python::vector_indexing_suite<std::vector<Boundary<Dim<2> >*> >())
    .def("__iter__", boost::python::iterator<std::vector<Boundary<Dim<2> >*>,
         return_internal_reference<> >())
    ;

  class_<std::vector<Boundary<Dim<3> >*> >("vector_of_Boundary3d_ptr")
    .def(boost::python::vector_indexing_suite<std::vector<Boundary<Dim<3> >*> >())
    .def("__iter__", boost::python::iterator<std::vector<Boundary<Dim<3> >*>,
         return_internal_reference<> >())
    ;

}

}
