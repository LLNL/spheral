#include <boost/python.hpp>
#include "Utilities/nodeOrdering.hh"
#include "Utilities/KeyTraits.hh"
#include "Geometry/Dimension.hh"
#include "Field/FieldList.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Expose the methods for computing the node ordering corresponding to sorting
// the given set of keys.
//------------------------------------------------------------------------------
void wrapNodeOrdering() {
  def("nodeOrdering", &nodeOrdering<Dim<1>, KeyTraits::Key>, "Compute the [1,n) numbering corresponding to sorting the given set of keys");
  def("nodeOrdering", &nodeOrdering<Dim<2>, KeyTraits::Key>, "Compute the [1,n) numbering corresponding to sorting the given set of keys");
  def("nodeOrdering", &nodeOrdering<Dim<3>, KeyTraits::Key>, "Compute the [1,n) numbering corresponding to sorting the given set of keys");
}

}
