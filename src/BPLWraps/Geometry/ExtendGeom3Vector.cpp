#include <string>
#include <sstream>
using namespace std;

#include "ExtendGeom3Vector.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Indexing
//------------------------------------------------------------------------------
double
index3Vector(Geom3Vector* self, int i) {
  return self->operator()(i);
}

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list
vector3Elements(Geom3Vector* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

//------------------------------------------------------------------------------
// Nice string printing.
//------------------------------------------------------------------------------
const char*
printGeom3Vector(Geom3Vector* self) {
  std::stringstream buffer;
  buffer << "3Vector("
         << self->x() << " "
         << self->y() << " "
         << self->z() << ")"
         << ends;
  return buffer.str().c_str();
}

}
