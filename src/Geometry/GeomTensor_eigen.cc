//---------------------------------Spheral++----------------------------------//
// GeomTensor -- the full tensor class.
//----------------------------------------------------------------------------//
#include "GeomTensor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
template<> const unsigned GeomTensor<1>::nDimensions = 1;
template<> const unsigned GeomTensor<1>::numElements = 1;
template<> const GeomTensor<1> GeomTensor<1>::zero = GeomTensor<1>(0.0);
template<> const GeomTensor<1> GeomTensor<1>::one = GeomTensor<1>(1.0);

template<> const unsigned GeomTensor<2>::nDimensions = 2;
template<> const unsigned GeomTensor<2>::numElements = 4;
template<> const GeomTensor<2> GeomTensor<2>::zero = GeomTensor<2>(0.0, 0.0,
                                                                   0.0, 0.0);
template<> const GeomTensor<2> GeomTensor<2>::one = GeomTensor<2>(1.0, 0.0,
                                                                  0.0, 1.0);

template<> const unsigned GeomTensor<3>::nDimensions = 3;
template<> const unsigned GeomTensor<3>::numElements = 9;
template<> const GeomTensor<3> GeomTensor<3>::zero = GeomTensor<3>(0.0, 0.0, 0.0,
                                                                   0.0, 0.0, 0.0,
                                                                   0.0, 0.0, 0.0);
template<> const GeomTensor<3> GeomTensor<3>::one = GeomTensor<3>(1.0, 0.0, 0.0,
                                                                  0.0, 1.0, 0.0,
                                                                  0.0, 0.0, 1.0);

}
