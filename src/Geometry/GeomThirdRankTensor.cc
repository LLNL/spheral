//---------------------------------Spheral++----------------------------------//
// GeomThirdRankTensor -- Geometric Tensor (rank 3) Class.
//
// This is a very simple, limited functionality rank 3 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim.
//
// Created by JMO, Thu Jul 17 17:05:58 PDT 2008
//----------------------------------------------------------------------------//
#include "GeomThirdRankTensor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
template<> const unsigned GeomThirdRankTensor<1>::numElements = 1;
template<> const GeomThirdRankTensor<1> GeomThirdRankTensor<1>::zero = GeomThirdRankTensor<1>(0.0);

template<> const unsigned GeomThirdRankTensor<2>::numElements = 8;
template<> const GeomThirdRankTensor<2> GeomThirdRankTensor<2>::zero = GeomThirdRankTensor<2>(0.0);

template<> const unsigned GeomThirdRankTensor<3>::numElements = 27;
template<> const GeomThirdRankTensor<3> GeomThirdRankTensor<3>::zero = GeomThirdRankTensor<3>(0.0);

}
