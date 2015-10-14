//---------------------------------Spheral++----------------------------------//
// invertRankNTensor -- Overloaded method to compute the inverse of our
// specialized rank N tensors.
//
// Created by JMO, Mon Oct 12 14:16:42 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_invertRankNTensor_hh__
#define __Spheral_invertRankNTensor_hh__

namespace Spheral {

// For now we only support a subset of the possible tensors.
template<typename Tensor>
Tensor invertRankNTensor(const Tensor& tensor);

}

#endif
