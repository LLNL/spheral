#include <boost/python.hpp>
#include "Utilities/initializeTau.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Wrap the initializeTau method.
//------------------------------------------------------------------------------
void wrapInitializeTau() {
  def("initializeTau", &initializeTau);
}

}
