#include "Geometry/Dimension.hh"

namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
WendlandC4Kernel< Dim<1> >::WendlandC4Kernel():
  Kernel<Dim<1>, WendlandC4Kernel< Dim<1> > >() {
  setVolumeNormalization(2.0/3.0);
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
WendlandC4Kernel< Dim<2> >::WendlandC4Kernel():
  Kernel<Dim<2>, WendlandC4Kernel< Dim<2> > >() {
  setVolumeNormalization(10.0/(7.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
WendlandC4Kernel< Dim<3> >::WendlandC4Kernel():
  Kernel<Dim<3>, WendlandC4Kernel< Dim<3> > >() {
  setVolumeNormalization(0.72973/M_PI);
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
WendlandC4Kernel<Dimension>::~WendlandC4Kernel() {
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC4Kernel<Dimension>::kernelValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 2.0) {
    double eta2 = etaMagnitude*etaMagnitude;
    return this->volumeNormalization()*Hdet*(pow(1.0-0.5*etaMagnitude,6.0)*
                                             (1.0+3.0*etaMagnitude+5.8333*eta2));
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
WendlandC4Kernel<Dimension>::gradValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 2.0) {
    return  this->volumeNormalization()*Hdet*(-3.0*pow(1.0-0.5*etaMagnitude,5.0)*
                                              (1.0+3.0*etaMagnitude+5.8333*eta2) +
                                              pow(1.0-0.5*etaMagnitude,6.0)*
                                              (3.0+2.91667*etaMagnitude));
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
WendlandC4Kernel<Dimension>::grad2Value(double etaMagnitude, double Hdet) const {
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
}
