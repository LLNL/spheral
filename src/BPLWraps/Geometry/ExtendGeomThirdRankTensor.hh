#include "Geometry/GeomThirdRankTensor.hh"
#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomVector.hh"

#include "BPLWraps/iteratorsAsList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Provide int, int, int indexing into tensors.
//------------------------------------------------------------------------------
double indexGeomThirdRankTensor1d(GeomThirdRankTensor<1>* self, int i, int j, int k);
double indexGeomThirdRankTensor2d(GeomThirdRankTensor<2>* self, int i, int j, int k);
double indexGeomThirdRankTensor3d(GeomThirdRankTensor<3>* self, int i, int j, int k);

void setGeomThirdRankTensor1d(GeomThirdRankTensor<1>* self, int i, int j, int k, double val);
void setGeomThirdRankTensor2d(GeomThirdRankTensor<2>* self, int i, int j, int k, double val);
void setGeomThirdRankTensor3d(GeomThirdRankTensor<3>* self, int i, int j, int k, double val);

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list thirdRankTensorElements1d(GeomThirdRankTensor<1>* self);
boost::python::list thirdRankTensorElements2d(GeomThirdRankTensor<2>* self);
boost::python::list thirdRankTensorElements3d(GeomThirdRankTensor<3>* self);

//------------------------------------------------------------------------------
// A nicer string representation of tensors.
//------------------------------------------------------------------------------
const char* printGeomThirdRankTensor1d(GeomThirdRankTensor<1>* self);
const char* printGeomThirdRankTensor2d(GeomThirdRankTensor<2>* self);
const char* printGeomThirdRankTensor3d(GeomThirdRankTensor<3>* self);

}
