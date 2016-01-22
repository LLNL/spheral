//---------------------------------Spheral++----------------------------------//
// GeomFourthRankTensor -- Geometric Tensor (rank 4) Class.
//
// This is a very simple, limited functionality rank 4 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim x nDim.
//
// Created by JMO, Sun Oct 11 22:10:34 PDT 2015
//----------------------------------------------------------------------------//
#include "GeomFourthRankTensor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
template<> const unsigned GeomFourthRankTensor<1>::numElements = 1;
template<> const GeomFourthRankTensor<1> GeomFourthRankTensor<1>::zero = GeomFourthRankTensor<1>(0.0);

template<> const unsigned GeomFourthRankTensor<2>::numElements = 16;
template<> const GeomFourthRankTensor<2> GeomFourthRankTensor<2>::zero = GeomFourthRankTensor<2>(0.0);

template<> const unsigned GeomFourthRankTensor<3>::numElements = 81;
template<> const GeomFourthRankTensor<3> GeomFourthRankTensor<3>::zero = GeomFourthRankTensor<3>(0.0);

}
