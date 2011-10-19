#include "Distributed/DomainNode.hh"
#include "Geometry/Dimension.hh"

#include "STLVectorWrap.hh"

// BPL includes
#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;
using namespace Spheral::PartitionSpace;

//------------------------------------------------------------------------------
// Wrap std::vector of Boundary.
//------------------------------------------------------------------------------
void
wrapSTLVectorDistributed() {

  class_<std::vector<DomainNode<Dim<2> > > >("vector_of_DomainNode2d")
    .def(boost::python::vector_indexing_suite<std::vector<DomainNode<Dim<2> > > >())
    ;

  class_<std::vector<DomainNode<Dim<3> > > >("vector_of_DomainNode3d")
    .def(boost::python::vector_indexing_suite<std::vector<DomainNode<Dim<3> > > >())
    ;

}

}
