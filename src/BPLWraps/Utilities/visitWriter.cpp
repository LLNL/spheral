#include <boost/python.hpp>
#include "ExtendVisitWriter.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Expose the routines for writing Visit mesh files.
//------------------------------------------------------------------------------
void wrapVisitWriter() {
  def("writeRectilinearMesh", &writeRectilinearMesh<Dim<2> >, "Write a 2-D Visit rectilinear mesh file from the sampled fields.");
  def("writeRectilinearMesh", &writeRectilinearMesh<Dim<3> >, "Write a 3-D Visit rectilinear mesh file from the sampled fields.");
}

}
