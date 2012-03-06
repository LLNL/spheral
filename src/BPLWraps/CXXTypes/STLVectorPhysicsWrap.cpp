#include "Physics/Physics.hh"
#include "Geometry/Dimension.hh"

// BPL includes
#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;
using namespace Spheral::PhysicsSpace;
using namespace Spheral::NodeSpace;

//------------------------------------------------------------------------------
// Wrap std::vector of Physics.
//------------------------------------------------------------------------------
void
wrapSTLVectorPhysics() {

  class_<std::vector<Physics<Dim<1> >*> >("vector_of_Physics1d_ptr")
    .def(boost::python::vector_indexing_suite<std::vector<Physics<Dim<1> >*> >())
    .def("__iter__", boost::python::iterator<std::vector<Physics<Dim<1> >*>,
         return_internal_reference<> >())
    ;

  class_<std::vector<Physics<Dim<2> >*> >("vector_of_Physics2d_ptr")
    .def(boost::python::vector_indexing_suite<std::vector<Physics<Dim<2> >*> >())
    .def("__iter__", boost::python::iterator<std::vector<Physics<Dim<2> >*>,
         return_internal_reference<> >())
    ;

  class_<std::vector<Physics<Dim<3> >*> >("vector_of_Physics3d_ptr")
    .def(boost::python::vector_indexing_suite<std::vector<Physics<Dim<3> >*> >())
    .def("__iter__", boost::python::iterator<std::vector<Physics<Dim<3> >*>,
         return_internal_reference<> >())
    ;

  // We also wrap the vector<pair<NodeList*, string>> here since it's used in the 
  // Physics state.  Strictly speaking this doesn't really belong in this function.
  class_<std::vector<std::pair<const NodeList<Dim<1> >*, std::string> > >("vector_of_pair_NodeList1d_string")
    .def(boost::python::vector_indexing_suite<std::vector<std::pair<NodeList<Dim<1> >*, std::string> > >())
    .def("__iter__", boost::python::iterator<std::vector<std::pair<const NodeList<Dim<1> >*, std::string> >,
         return_internal_reference<> >())
    ;

  class_<std::vector<std::pair<const NodeList<Dim<2> >*, std::string> > >("vector_of_pair_NodeList2d_string")
    .def(boost::python::vector_indexing_suite<std::vector<std::pair<NodeList<Dim<2> >*, std::string> > >())
    .def("__iter__", boost::python::iterator<std::vector<std::pair<const NodeList<Dim<2> >*, std::string> >,
         return_internal_reference<> >())
    ;

  class_<std::vector<std::pair<const NodeList<Dim<3> >*, std::string> > >("vector_of_pair_NodeList3d_string")
    .def(boost::python::vector_indexing_suite<std::vector<std::pair<NodeList<Dim<3> >*, std::string> > >())
    .def("__iter__", boost::python::iterator<std::vector<std::pair<const NodeList<Dim<3> >*, std::string> >,
         return_internal_reference<> >())
    ;
}

}
