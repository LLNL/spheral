#include "QuinticSplineKernel.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {
namespace KernelSpace {

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
QuinticSplineKernel<Dimension>::~QuinticSplineKernel() {
}

//------------------------------------------------------------------------------
// Return the kernel weight for a given normalized distance.
//------------------------------------------------------------------------------
template<typename Dimension>
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
