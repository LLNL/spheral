//---------------------------------Spheral++----------------------------------//
// GeomSymmetricTensor -- the symmetric tensor class
//----------------------------------------------------------------------------//
#include <cmath>
#include <limits>
#include <float.h>
#include <vector>

#include "GeomSymmetricTensor.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class GeomSymmetricTensor<1>;
template class GeomSymmetricTensor<2>;
template class GeomSymmetricTensor<3>;

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
template<> const unsigned GeomSymmetricTensor<1>::nDimensions = 1;
template<> const unsigned GeomSymmetricTensor<1>::numElements = 1;
template<> const GeomSymmetricTensor<1> GeomSymmetricTensor<1>::zero = GeomSymmetricTensor<1>(0.0);
template<> const GeomSymmetricTensor<1> GeomSymmetricTensor<1>::one = GeomSymmetricTensor<1>(1.0);

template<> const unsigned GeomSymmetricTensor<2>::nDimensions = 2;
template<> const unsigned GeomSymmetricTensor<2>::numElements = 4;
template<> const GeomSymmetricTensor<2> GeomSymmetricTensor<2>::zero = GeomSymmetricTensor<2>(0.0, 0.0,
                                                                                              0.0, 0.0);
template<> const GeomSymmetricTensor<2> GeomSymmetricTensor<2>::one = GeomSymmetricTensor<2>(1.0, 0.0,
                                                                                             0.0, 1.0);

template<> const unsigned GeomSymmetricTensor<3>::nDimensions = 3;
template<> const unsigned GeomSymmetricTensor<3>::numElements = 9;
template<> const GeomSymmetricTensor<3> GeomSymmetricTensor<3>::zero = GeomSymmetricTensor<3>(0.0, 0.0, 0.0,
                                                                                              0.0, 0.0, 0.0,
                                                                                              0.0, 0.0, 0.0);
template<> const GeomSymmetricTensor<3> GeomSymmetricTensor<3>::one = GeomSymmetricTensor<3>(1.0, 0.0, 0.0,
                                                                                             0.0, 1.0, 0.0,
                                                                                             0.0, 0.0, 1.0);

template<> const double GeomSymmetricTensor<1>::onethird = 1.0/3.0;
template<> const double GeomSymmetricTensor<2>::onethird = 1.0/3.0;
template<> const double GeomSymmetricTensor<3>::onethird = 1.0/3.0;

template<> const double GeomSymmetricTensor<1>::sqrt3 = std::sqrt(3.0);
template<> const double GeomSymmetricTensor<2>::sqrt3 = std::sqrt(3.0);
template<> const double GeomSymmetricTensor<3>::sqrt3 = std::sqrt(3.0);

}

