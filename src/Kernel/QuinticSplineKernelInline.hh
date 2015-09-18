#include "Geometry/Dimension.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
template<>
inline
QuinticSplineKernel< Dim<1> >::QuinticSplineKernel():
  Kernel<Dim<1>, QuinticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(3.0/48.0);
  setKernelExtent(2.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}

template<>
inline
QuinticSplineKernel< Dim<2> >::QuinticSplineKernel():
  Kernel<Dim<2>, QuinticSplineKernel< Dim<2> > >() {
  setVolumeNormalization(3.0/(16.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}

template<>
inline
QuinticSplineKernel< Dim<3> >::QuinticSplineKernel():
  Kernel<Dim<3>, QuinticSplineKernel< Dim<3> > >() {
  setVolumeNormalization(7.0/(40.0*M_PI));
  setKernelExtent(2.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
QuinticSplineKernel<Dimension>::~QuinticSplineKernel() {
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
QuinticSplineKernel<Dimension>::kernelValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaMagnitude < 1.0) {
    return this->volumeNormalization()*Hdet*(FastMath::pow5(2.0 - etaMagnitude) - 
                                             16.0*FastMath::pow5(1.0 - etaMagnitude));
  } else if (etaMagnitude < 2.0) {
    return this->volumeNormalization()*Hdet*(FastMath::pow5(2.0 - etaMagnitude));
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
QuinticSplineKernel<Dimension>::gradValue(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaMagnitude < 1.0) {
    return this->volumeNormalization()*Hdet*(-5.0*FastMath::pow4(2.0 - etaMagnitude) +
                                             80.0*FastMath::pow4(1.0 - etaMagnitude));
  } else if (etaMagnitude < 2.0) {
    return this->volumeNormalization()*Hdet*(-5.0*FastMath::pow4(2.0 - etaMagnitude));
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
QuinticSplineKernel<Dimension>::grad2Value(double etaMagnitude, double Hdet) const {
  REQUIRE(etaMagnitude >= 0.0);
  REQUIRE(Hdet >= 0.0);
  if (etaMagnitude < 1.0) {
    return this->volumeNormalization()*Hdet*(20.0*FastMath::pow3(2.0 - etaMagnitude) - 
                                             320.0*FastMath::pow3(1.0 - etaMagnitude));
  } else if (etaMagnitude < 2.0) {
    return this->volumeNormalization()*Hdet*(20.0*FastMath::pow3(2.0 - etaMagnitude));
  } else {
    return 0.0;
  }
}

}
}
