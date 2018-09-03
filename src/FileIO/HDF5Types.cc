#include "HDF5Types.hh"
#include "Geometry/Dimension.hh"

#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// A function to define the static compound data types for HDF5 used by Spheral.
//------------------------------------------------------------------------------
int initializeSpheralH5Types() {

  // Initialize the storage for each of the Vector types.
  H5Vector1d = CompType(sizeof(Dim<1>::Vector));
  H5Vector2d = CompType(sizeof(Dim<2>::Vector));
  H5Vector3d = CompType(sizeof(Dim<3>::Vector));

  // 1-D Vector.
  H5Vector1d.insertMember("x", 0*sizeof(double), PredType::NATIVE_DOUBLE);

  // 2-D Vector.
  H5Vector2d.insertMember("x", 0*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Vector2d.insertMember("y", 1*sizeof(double), PredType::NATIVE_DOUBLE);

  // 3-D Vector.
  H5Vector3d.insertMember("x", 0*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Vector3d.insertMember("y", 1*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Vector3d.insertMember("z", 2*sizeof(double), PredType::NATIVE_DOUBLE);

  // Initialize the storage for the Tensor types.
  H5Tensor1d = CompType(sizeof(Dim<1>::Tensor));
  H5Tensor2d = CompType(sizeof(Dim<2>::Tensor));
  H5Tensor3d = CompType(sizeof(Dim<3>::Tensor));

  // 1-D Tensor.
  H5Tensor1d.insertMember("xx", 0*sizeof(double), PredType::NATIVE_DOUBLE);

  // 2-D Tensor.
  H5Tensor2d.insertMember("xx", 0*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor2d.insertMember("xy", 1*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor2d.insertMember("yx", 2*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor2d.insertMember("yy", 3*sizeof(double), PredType::NATIVE_DOUBLE);

  // 3-D Tensor.
  H5Tensor3d.insertMember("xx", 0*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("xy", 1*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("xz", 2*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("yx", 3*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("yy", 4*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("yz", 5*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("zx", 6*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("zy", 7*sizeof(double), PredType::NATIVE_DOUBLE);
  H5Tensor3d.insertMember("zz", 8*sizeof(double), PredType::NATIVE_DOUBLE);

  // Initialize the storage for the Symmetric Tensor types.
  H5SymTensor1d = CompType(sizeof(Dim<1>::SymTensor));
  H5SymTensor2d = CompType(sizeof(Dim<2>::SymTensor));
  H5SymTensor3d = CompType(sizeof(Dim<3>::SymTensor));

  // 1-D SymTensor.
  H5SymTensor1d.insertMember("xx", 0*sizeof(double), PredType::NATIVE_DOUBLE);

  // 2-D SymTensor.
  H5SymTensor2d.insertMember("xx", 0*sizeof(double), PredType::NATIVE_DOUBLE);
  H5SymTensor2d.insertMember("xy", 1*sizeof(double), PredType::NATIVE_DOUBLE);
  H5SymTensor2d.insertMember("yy", 2*sizeof(double), PredType::NATIVE_DOUBLE);

  // 3-D SymTensor.
  H5SymTensor3d.insertMember("xx", 0*sizeof(double), PredType::NATIVE_DOUBLE);
  H5SymTensor3d.insertMember("xy", 1*sizeof(double), PredType::NATIVE_DOUBLE);
  H5SymTensor3d.insertMember("xz", 2*sizeof(double), PredType::NATIVE_DOUBLE);
  H5SymTensor3d.insertMember("yy", 3*sizeof(double), PredType::NATIVE_DOUBLE);
  H5SymTensor3d.insertMember("yz", 4*sizeof(double), PredType::NATIVE_DOUBLE);
  H5SymTensor3d.insertMember("zz", 5*sizeof(double), PredType::NATIVE_DOUBLE);

}

//------------------------------------------------------------------------------
// Call the function to fill in the static H5 Spheral types.
//------------------------------------------------------------------------------
const int initTypes = initializeSpheralH5Types();

}
