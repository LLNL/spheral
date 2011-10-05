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
// Explicit instantiation.
//------------------------------------------------------------------------------
template class GeomThirdRankTensor<1>;
template class GeomThirdRankTensor<2>;
template class GeomThirdRankTensor<3>;

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
template<> const unsigned GeomThirdRankTensor<1>::nDimensions = 1;
template<> const unsigned GeomThirdRankTensor<1>::numElements = 1;
template<> const unsigned GeomThirdRankTensor<1>::nDim2       = 1;
template<> const GeomThirdRankTensor<1> GeomThirdRankTensor<1>::zero = GeomThirdRankTensor<1>(0.0);

template<> const unsigned GeomThirdRankTensor<2>::nDimensions = 2;
template<> const unsigned GeomThirdRankTensor<2>::numElements = 8;
template<> const unsigned GeomThirdRankTensor<2>::nDim2       = 4;
template<> const GeomThirdRankTensor<2> GeomThirdRankTensor<2>::zero = GeomThirdRankTensor<2>(0.0);

template<> const unsigned GeomThirdRankTensor<3>::nDimensions = 3;
template<> const unsigned GeomThirdRankTensor<3>::numElements = 27;
template<> const unsigned GeomThirdRankTensor<3>::nDim2       = 9;
template<> const GeomThirdRankTensor<3> GeomThirdRankTensor<3>::zero = GeomThirdRankTensor<3>(0.0);

}
