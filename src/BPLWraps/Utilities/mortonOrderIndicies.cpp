#include <boost/python.hpp>
#include "Utilities/mortonOrderIndicies.hh"
#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Expose the routines for computing the morton ordering indicies.
//------------------------------------------------------------------------------
void wrapMortonOrderIndicies() {
  def("mortonOrderIndicies", &mortonOrderIndicies<Dim<1> >, "Compute the morton ordering indicies for all nodes in the DataBase.");
  def("mortonOrderIndicies", &mortonOrderIndicies<Dim<2> >, "Compute the morton ordering indicies for all nodes in the DataBase.");
  def("mortonOrderIndicies", &mortonOrderIndicies<Dim<3> >, "Compute the morton ordering indicies for all nodes in the DataBase.");
}

}
