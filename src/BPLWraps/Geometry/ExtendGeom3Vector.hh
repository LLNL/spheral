#include "Geometry/Geom3Vector.hh"
#include "BPLWraps/iteratorsAsList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Indexing
//------------------------------------------------------------------------------
double
index3Vector(Geom3Vector* self, int i);

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list
vector3Elements(Geom3Vector* self);

//------------------------------------------------------------------------------
// Nice string printing.
//------------------------------------------------------------------------------
const char*
printGeom3Vector(Geom3Vector* self);

}
