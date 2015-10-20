//---------------------------------Spheral++----------------------------------//
// GeomTensor -- the full tensor class.
//----------------------------------------------------------------------------//
#include <cmath>
#include <limits.h>
#include <float.h>

#include "GeomTensor.hh"
#include "findEigenValues3.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/rotationMatrix.hh"

using namespace std;

namespace Spheral {

//------------------------------------------------------------------------------
// Find the eigenvalues of a tensor.
//------------------------------------------------------------------------------
template<>
GeomVector<1>
GeomTensor<1>::eigenValues() const {
  return GeomVector<1>(this->mxx);
}

//----------------------------------------------------------------------
template<>
GeomVector<2>
GeomTensor<2>::eigenValues() const {
  const double b = Trace();
  const double c = Determinant();
  const double q = 0.5*(b + sgn(b)*sqrt(std::max(0.0, b*b - 4.0*c))) + 1.0e-10*sgn(b);
  CHECK(q != 0.0);
  return GeomVector<2>(q, c/q);
}

//----------------------------------------------------------------------
template<>
GeomVector<3>
GeomTensor<3>::eigenValues() const {
  return findEigenValues3<GeomTensor<3> >(*this);
}

// //------------------------------------------------------------------------------
// // Explicit instantiation.
// //------------------------------------------------------------------------------
// template class GeomTensor<1>;
// template class GeomTensor<2>;
// template class GeomTensor<3>;

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
