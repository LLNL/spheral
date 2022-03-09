#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
BSplineKernel< Dim<1> >::BSplineKernel():
  Kernel<Dim<1>, BSplineKernel< Dim<1> > >() {
  setVolumeNormalization(2.0/3.0);
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
BSplineKernel< Dim<2> >::BSplineKernel():
  Kernel<Dim<2>, BSplineKernel< Dim<2> > >() {
  setVolumeNormalization(10.0/(7.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
BSplineKernel< Dim<3> >::BSplineKernel():
  Kernel<Dim<3>, BSplineKernel< Dim<3> > >() {
  setVolumeNormalization(1.0/M_PI);
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
BSplineKernel<Dimension>::~BSplineKernel() {
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
BSplineKernel<Dimension>::kernelValue(double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 1.0) {
    double eta2 = etaMagnitude*etaMagnitude;
    return this->volumeNormalization()*Hdet*(1.0 - 1.5*eta2 + 0.75*eta2*etaMagnitude);
  } else if (etaMagnitude < 2.0) {
    return this->volumeNormalization()*Hdet*0.25*FastMath::pow3(2.0 - etaMagnitude);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the gradient value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
BSplineKernel<Dimension>::gradValue(double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 1.0) {
    return  -this->volumeNormalization()*Hdet*(3.0 - 2.25*etaMagnitude)*etaMagnitude;
  } else if (etaMagnitude < 2.0) {
    return -this->volumeNormalization()*Hdet*0.75*FastMath::pow2(2.0 - etaMagnitude);
  } else {
    return 0.0;
  }
}

//------------------------------------------------------------------------------
// Return the second derivative value for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
BSplineKernel<Dimension>::grad2Value(double etaMagnitude, const double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 1.0) {
    return  -this->volumeNormalization()*Hdet*(3 - 4.5*etaMagnitude);
  } else if (etaMagnitude < 2.0) {
    return this->volumeNormalization()*Hdet*1.5*(2 - etaMagnitude);
  } else {
    return 0.0;
  }
}

}
