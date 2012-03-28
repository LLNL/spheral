#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomVector.hh"
#include "Geometry/EigenStruct.hh"

#include "BPLWraps/iteratorsAsList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Provide int, int indexing into tensors.
//------------------------------------------------------------------------------
double
indexTensor1d(GeomTensor<1>* self, int row, int col);

double
indexSymTensor1d(GeomSymmetricTensor<1>* self, int row, int col);

double
indexTensor2d(GeomTensor<2>* self, int row, int col);

double
indexSymTensor2d(GeomSymmetricTensor<2>* self, int row, int col);

double
indexTensor3d(GeomTensor<3>* self, int row, int col);

double
indexSymTensor3d(GeomSymmetricTensor<3>* self, int row, int col);

//------------------------------------------------------------------------------
// List access to elements.
//------------------------------------------------------------------------------
boost::python::list
tensorElements1d(GeomTensor<1>* self);

boost::python::list
tensorElements2d(GeomTensor<2>* self);

boost::python::list
tensorElements3d(GeomTensor<3>* self);

boost::python::list
symTensorElements1d(GeomSymmetricTensor<1>* self);

boost::python::list
symTensorElements2d(GeomSymmetricTensor<2>* self);

boost::python::list
symTensorElements3d(GeomSymmetricTensor<3>* self);

//------------------------------------------------------------------------------
// A nicer string representation of tensors.
//------------------------------------------------------------------------------
const char*
printTensor1d(GeomTensor<1>* self);

const char*
printTensor2d(GeomTensor<2>* self);

const char*
printTensor3d(GeomTensor<3>* self);

const char*
printSymTensor1d(GeomSymmetricTensor<1>* self);

const char*
printSymTensor2d(GeomSymmetricTensor<2>* self);

const char*
printSymTensor3d(GeomSymmetricTensor<3>* self);

}
