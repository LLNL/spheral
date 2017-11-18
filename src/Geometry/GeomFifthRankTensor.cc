//---------------------------------Spheral++----------------------------------//
// GeomFifthRankTensor -- Geometric Tensor (rank 4) Class.
//
// This is a very simple, limited functionality rank 4 tensor.  Assumes
// the thing is dimensioned nDim x nDim x nDim x nDim.
//
// Created by JMO, Tue Oct 13 15:24:56 PDT 2015
//----------------------------------------------------------------------------//
#include "GeomFifthRankTensor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
template<> const unsigned GeomFifthRankTensor<1>::numElements = 1;
template<> const GeomFifthRankTensor<1> GeomFifthRankTensor<1>::zero = GeomFifthRankTensor<1>(0.0);

template<> const unsigned GeomFifthRankTensor<2>::numElements = 32;
template<> const GeomFifthRankTensor<2> GeomFifthRankTensor<2>::zero = GeomFifthRankTensor<2>(0.0);

template<> const unsigned GeomFifthRankTensor<3>::numElements = 243;
template<> const GeomFifthRankTensor<3> GeomFifthRankTensor<3>::zero = GeomFifthRankTensor<3>(0.0);

}
