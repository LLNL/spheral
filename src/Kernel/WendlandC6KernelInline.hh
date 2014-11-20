#include "Geometry/Dimension.hh"

namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
WendlandC6Kernel< Dim<1> >::WendlandC6Kernel():
  Kernel<Dim<1>, WendlandC6Kernel< Dim<1> > >() {
  setVolumeNormalization(99.0/74.0);
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
WendlandC6Kernel< Dim<2> >::WendlandC6Kernel():
  Kernel<Dim<2>, WendlandC6Kernel< Dim<2> > >() {
  setVolumeNormalization(6435.0/(3856.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
WendlandC6Kernel< Dim<3> >::WendlandC6Kernel():
  Kernel<Dim<3>, WendlandC6Kernel< Dim<3> > >() {
  setVolumeNormalization(45045.0/(31456.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
WendlandC6Kernel<Dimension>::~WendlandC6Kernel() {
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
WendlandC6Kernel<Dimension>::kernelValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 2.0) {
    double eta2 = etaMagnitude*etaMagnitude;
    return this->volumeNormalization()*Hdet*(pow(1.0-0.5*etaMagnitude,8)*
                                             (1.0+4.0*etaMagnitude+(25.0/2.0)*eta2+18.0*etaMagnitude*eta2));
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
WendlandC6Kernel<Dimension>::gradValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 2.0) {
    double eta2 = etaMagnitude*etaMagnitude;
    return  this->volumeNormalization()*Hdet*((1.0/256.0)*pow(etaMagnitude-2.0,7)*
                                              etaMagnitude*(198.0*eta2+17*etaMagnitude-14.0));
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
WendlandC6Kernel<Dimension>::grad2Value(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);

  if (etaMagnitude < 1.0) {
    const double eta2 = etaMagnitude*etaMagnitude;
    return this->volumeNormalization()*Hdet*((1.0/256.0)*pow(etaMagnitude-2.0,6)*
                                               (1980.0*eta2*etaMagnitude-1035.0*eta2-180.0*etaMagnitude+28.0));
  } else {
    return 0.0;
  }
}

}
}
