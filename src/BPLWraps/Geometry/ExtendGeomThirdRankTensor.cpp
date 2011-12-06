#include <sstream>

#include "Geometry/GeomThirdRankTensor.hh"
#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomVector.hh"

#include "BPLWraps/iteratorsAsList.hh"

using namespace std;

namespace Spheral {

//------------------------------------------------------------------------------
// Provide int, int, int indexing into tensors. (get)
//------------------------------------------------------------------------------
double
indexGeomThirdRankTensor1d(GeomThirdRankTensor<1>* self, int i, int j, int k) {
  return self->operator()(i, j, k);
}

double
indexGeomThirdRankTensor2d(GeomThirdRankTensor<2>* self, int i, int j, int k) {
  return self->operator()(i, j, k);
}

double
indexGeomThirdRankTensor3d(GeomThirdRankTensor<3>* self, int i, int j, int k) {
  return self->operator()(i, j, k);
}

//------------------------------------------------------------------------------
// Provide int, int, int indexing into tensors. (set)
//------------------------------------------------------------------------------
void
setGeomThirdRankTensor1d(GeomThirdRankTensor<1>* self, int i, int j, int k, double val) {
  self->operator()(i, j, k) = val;
}

void
setGeomThirdRankTensor2d(GeomThirdRankTensor<2>* self, int i, int j, int k, double val) {
  self->operator()(i, j, k) = val;
}

void
setGeomThirdRankTensor3d(GeomThirdRankTensor<3>* self, int i, int j, int k, double val) {
  self->operator()(i, j, k) = val;
}

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list
thirdRankTensorElements1d(GeomThirdRankTensor<1>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
thirdRankTensorElements2d(GeomThirdRankTensor<2>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

boost::python::list
thirdRankTensorElements3d(GeomThirdRankTensor<3>* self) {
  return iteratorsAsListByValue(self->begin(),
                                self->end());
}

//------------------------------------------------------------------------------
// A nicer string representation of tensors.
//------------------------------------------------------------------------------
template<typename T>
inline
const char*
printGeomThirdRankTensor(T* self) {
  std::stringstream buffer;
  buffer << "ThirdRankTensor" << *self << ends;
  return buffer.str().c_str();
}

const char*
printGeomThirdRankTensor1d(GeomThirdRankTensor<1>* self) {
  return printGeomThirdRankTensor< GeomThirdRankTensor<1> >(self);
}

const char*
printGeomThirdRankTensor2d(GeomThirdRankTensor<2>* self) {
  return printGeomThirdRankTensor< GeomThirdRankTensor<2> >(self);
}

const char*
printGeomThirdRankTensor3d(GeomThirdRankTensor<3>* self) {
  return printGeomThirdRankTensor< GeomThirdRankTensor<3> >(self);
}

}
