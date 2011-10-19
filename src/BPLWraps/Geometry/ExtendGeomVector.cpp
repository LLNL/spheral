#include <string>
#include <sstream>
using namespace std;

#include "ExtendGeomVector.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Indexing
//------------------------------------------------------------------------------
double
indexVector1d(GeomVector<1>* self, int i) {
  return self->operator()(i);
}

double
indexVector2d(GeomVector<2>* self, int i) {
  return self->operator()(i);
}

double
indexVector3d(GeomVector<3>* self, int i) {
  return self->operator()(i);
}

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list
vectorElements1d(GeomVector<1>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
vectorElements2d(GeomVector<2>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
vectorElements3d(GeomVector<3>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

//------------------------------------------------------------------------------
// Nice string printing.
//------------------------------------------------------------------------------
template<int nDim>
inline
const char*
printGeomVector(GeomVector<nDim>* self) {
  std::stringstream buffer;
  buffer << "Vector(";
  for (int i = 0; i < nDim; ++i) {
    buffer << (*self)(i);
    if (i < nDim - 1) buffer << " ";
  }
  buffer << ")" << ends;
  return buffer.str().c_str();
}

const char*
printGeomVector1d(GeomVector<1>* self) {
  return printGeomVector<1>(self);
}

const char*
printGeomVector2d(GeomVector<2>* self) {
  return printGeomVector<2>(self);
}

const char*
printGeomVector3d(GeomVector<3>* self) {
  return printGeomVector<3>(self);
}

}
