#include <boost/python.hpp>
#include "Utilities/peanoHilbertOrderIndicies.hh"
#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Expose the routines for computing the Peano-Hilbert ordering indicies.
//------------------------------------------------------------------------------
void wrapPeanoHilbertOrderIndicies() {
  def("peanoHilbertOrderIndicies", &peanoHilbertOrderIndicies<Dim<1> >, "Compute the Peano-Hilbert ordering indicies for all nodes in the DataBase.");
  def("peanoHilbertOrderIndicies", &peanoHilbertOrderIndicies<Dim<2> >, "Compute the Peano-Hilbert ordering indicies for all nodes in the DataBase.");
  def("peanoHilbertOrderIndicies", &peanoHilbertOrderIndicies<Dim<3> >, "Compute the Peano-Hilbert ordering indicies for all nodes in the DataBase.");
}

}
