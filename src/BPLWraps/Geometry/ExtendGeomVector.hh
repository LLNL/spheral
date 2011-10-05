#include "Geometry/GeomVector.hh"
#include "BPLWraps/iteratorsAsList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Indexing
//------------------------------------------------------------------------------
double
indexVector1d(GeomVector<1>* self, int i);

double
indexVector2d(GeomVector<2>* self, int i);

double
indexVector3d(GeomVector<3>* self, int i);

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list
vectorElements1d(GeomVector<1>* self);

boost::python::list
vectorElements2d(GeomVector<2>* self);

boost::python::list
vectorElements3d(GeomVector<3>* self);

//------------------------------------------------------------------------------
// Nice string printing.
//------------------------------------------------------------------------------
const char*
printGeomVector1d(GeomVector<1>* self);

const char*
printGeomVector2d(GeomVector<2>* self);

const char*
printGeomVector3d(GeomVector<3>* self);

}
