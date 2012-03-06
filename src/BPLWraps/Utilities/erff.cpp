#include <boost/python.hpp>
#include "Utilities/erff.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Wrap the erff method.
//------------------------------------------------------------------------------
void wrapErff() {
  def("erff", &erff);
}

}
