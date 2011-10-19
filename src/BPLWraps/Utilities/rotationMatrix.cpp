#include <boost/python.hpp>
#include "Utilities/rotationMatrix.hh"
#include "Geometry/Dimension.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Expose the rotation matrix functions.
//------------------------------------------------------------------------------
void wrapRotationMatrix() {
  Dim<1>::Tensor (*rm1d)(const Dim<1>::Vector&) = &rotationMatrix;
  Dim<2>::Tensor (*rm2d)(const Dim<2>::Vector&) = &rotationMatrix;
  Dim<3>::Tensor (*rm3d)(const Dim<3>::Vector&) = &rotationMatrix;

  def("rotationMatrix", rm1d, "Compute a rotation matrix to transform into the frame with the x-axis aligned with the given unit vector.");
  def("rotationMatrix", rm2d, "Compute a rotation matrix to transform into the frame with the x-axis aligned with the given unit vector.");
  def("rotationMatrix", rm3d, "Compute a rotation matrix to transform into the frame with the x-axis aligned with the given unit vector.");
}

}
