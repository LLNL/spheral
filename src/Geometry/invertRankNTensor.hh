//---------------------------------Spheral++----------------------------------//
// invertRankNTensor -- Overloaded method to compute the inverse of our
// specialized rank N tensors.
//
// Created by JMO, Mon Oct 12 14:16:42 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_invertRankNTensor_hh__
#define __Spheral_invertRankNTensor_hh__

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

// General method.
template<typename Tensor>
Tensor invertRankNTensor(const Tensor& tensor);

// We might as well take advantage of our specialized lower-rank cases.
template<> inline Dim<1>::Tensor    invertRankNTensor(const Dim<1>::Tensor&    tensor) { return tensor.Inverse(); }
template<> inline Dim<1>::SymTensor invertRankNTensor(const Dim<1>::SymTensor& tensor) { return tensor.Inverse(); }
template<> inline Dim<2>::Tensor    invertRankNTensor(const Dim<2>::Tensor&    tensor) { return tensor.Inverse(); }
template<> inline Dim<2>::SymTensor invertRankNTensor(const Dim<2>::SymTensor& tensor) { return tensor.Inverse(); }
template<> inline Dim<3>::Tensor    invertRankNTensor(const Dim<3>::Tensor&    tensor) { return tensor.Inverse(); }
template<> inline Dim<3>::SymTensor invertRankNTensor(const Dim<3>::SymTensor& tensor) { return tensor.Inverse(); }

}

#endif
