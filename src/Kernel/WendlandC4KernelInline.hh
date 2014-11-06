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
  setVolumeNormalization(54.0/37.0);
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
WendlandC4Kernel< Dim<2> >::WendlandC4Kernel():
  Kernel<Dim<2>, WendlandC4Kernel< Dim<2> > >() {
  setVolumeNormalization(9.0/(5.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(1.0/3.0);
}

template<>
inline
WendlandC4Kernel< Dim<3> >::WendlandC4Kernel():
  Kernel<Dim<3>, WendlandC4Kernel< Dim<3> > >() {
  setVolumeNormalization(165.0/(112.0*M_PI));
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
    return this->volumeNormalization()*Hdet*(pow(1.0-0.5*etaMagnitude,6)*
                                             (1.0+3.0*etaMagnitude+(35.0/6.0)*eta2));
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
    return  this->volumeNormalization()*Hdet*((7.0/192.0)*pow(etaMagnitude-2.0,5)*
                                              etaMagnitude*(20.0*etaMagnitude-1.0));
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
    const double eta2 = etaMagnitude*etaMagnitude;
    return this->volumeNormalization()*Hdet*((7.0/96.0)*pow(etaMagnitude-2.0,4)*
                                               (70.0*eta2-43.0*etaMagnitude+1.0));
  } else {
    return 0.0;
  }
}

}
}
