#include "HDF5Traits.hh"
#include "HDF5Types.hh"
#include "Geometry/Dimension.hh"

#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// Set the H5Traits translation types.
//------------------------------------------------------------------------------
const CompType HDF5Traits<Dim<1>::Vector>::Type = H5Vector1d;
const CompType HDF5Traits<Dim<2>::Vector>::Type = H5Vector2d;
const CompType HDF5Traits<Dim<3>::Vector>::Type = H5Vector3d;
                                                  
const CompType HDF5Traits<Dim<1>::Tensor>::Type = H5Tensor1d;
const CompType HDF5Traits<Dim<2>::Tensor>::Type = H5Tensor2d;
const CompType HDF5Traits<Dim<3>::Tensor>::Type = H5Tensor3d;

const CompType HDF5Traits<Dim<1>::SymTensor>::Type = H5SymTensor1d;
const CompType HDF5Traits<Dim<2>::SymTensor>::Type = H5SymTensor2d;
const CompType HDF5Traits<Dim<3>::SymTensor>::Type = H5SymTensor3d;

}
