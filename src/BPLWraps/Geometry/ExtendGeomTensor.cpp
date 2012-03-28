#include <sstream>

#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomVector.hh"
#include "Geometry/EigenStruct.hh"

#include "BPLWraps/iteratorsAsList.hh"

using namespace std;

namespace Spheral {

//------------------------------------------------------------------------------
// Provide int, int indexing into tensors.
//------------------------------------------------------------------------------
double
indexTensor1d(GeomTensor<1>* self, int row, int col) {
  return self->operator()(row, col);
}

double
indexSymTensor1d(GeomSymmetricTensor<1>* self, int row, int col) {
  return self->operator()(row, col);
}

double
indexTensor2d(GeomTensor<2>* self, int row, int col) {
  return self->operator()(row, col);
}

double
indexSymTensor2d(GeomSymmetricTensor<2>* self, int row, int col) {
  return self->operator()(row, col);
}

double
indexTensor3d(GeomTensor<3>* self, int row, int col) {
  return self->operator()(row, col);
}

double
indexSymTensor3d(GeomSymmetricTensor<3>* self, int row, int col) {
  return self->operator()(row, col);
}

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list
tensorElements1d(GeomTensor<1>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
tensorElements2d(GeomTensor<2>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
tensorElements3d(GeomTensor<3>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
symTensorElements1d(GeomSymmetricTensor<1>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
symTensorElements2d(GeomSymmetricTensor<2>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
symTensorElements3d(GeomSymmetricTensor<3>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

//------------------------------------------------------------------------------
// A nicer string representation of tensors.
//------------------------------------------------------------------------------
template<typename TensorType>
inline
const char*
printGeomTensor(TensorType* self) {
  std::stringstream buffer;
  buffer << "Tensor(";
  for (int i = 0; i < self->nDimensions; ++i) {
    buffer << "(";
    for (int j = 0; j < self->nDimensions; ++j) {
      buffer << (*self)(i,j);
      if (j < self->nDimensions - 1) buffer << " ";
    }
    buffer << ")";
  }
  buffer << ")" << ends;
  return buffer.str().c_str();
}

const char*
printTensor1d(GeomTensor<1>* self) {
  return printGeomTensor< GeomTensor<1> >(self);
}

const char*
printTensor2d(GeomTensor<2>* self) {
  return printGeomTensor< GeomTensor<2> >(self);
}

const char*
printTensor3d(GeomTensor<3>* self) {
  return printGeomTensor< GeomTensor<3> >(self);
}

const char*
printSymTensor1d(GeomSymmetricTensor<1>* self) {
  return printGeomTensor< GeomSymmetricTensor<1> >(self);
}

const char*
printSymTensor2d(GeomSymmetricTensor<2>* self) {
  return printGeomTensor< GeomSymmetricTensor<2> >(self);
}

const char*
printSymTensor3d(GeomSymmetricTensor<3>* self) {
  return printGeomTensor< GeomSymmetricTensor<3> >(self);
}

}
