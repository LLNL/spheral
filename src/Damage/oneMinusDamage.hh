#ifndef __Spheral__oneMinusDamage__
#define __Spheral__oneMinusDamage__

#include "Geometry/Dimension.hh"

namespace Spheral {

inline
Dim<1>::SymTensor
oneMinusDamage(const Dim<1>::SymTensor& Di) {
  REQUIRE(Di.xx() >= 0.0);
  Dim<1>::SymTensor result(std::max(0.0, std::min(1.0, 1.0 - Di.xx())));
  ENSURE(result.xx() >= 0.0 &&
         result.xx() <= 1.0);
  return result;
}

inline
Dim<2>::SymTensor
oneMinusDamage(const Dim<2>::SymTensor& Di) {
  const Dim<2>::SymTensor::EigenStructType eigen = Di.eigenVectors();
  REQUIRE(eigen.eigenValues.minElement() >= -1e-5);
  Dim<2>::SymTensor result(std::max(1.0e-10, std::min(1.0, 1.0 - eigen.eigenValues.x())),
                           0.0,
                           0.0,
                           std::max(1.0e-10, std::min(1.0, 1.0 - eigen.eigenValues.y())));
  result.rotationalTransform(eigen.eigenVectors);
  ENSURE(result.eigenValues().minElement() >= -1e-5 &&
         result.eigenValues().maxElement() <= 1.0 + 1e-5);
  return result;
}

inline
Dim<3>::SymTensor
oneMinusDamage(const Dim<3>::SymTensor& Di) {
  const Dim<3>::SymTensor::EigenStructType eigen = Di.eigenVectors();
  REQUIRE(fuzzyGreaterThanOrEqual(eigen.eigenValues.minElement(), 0.0, 1.0e-8));
  Dim<3>::SymTensor result(std::max(1.0e-4, std::min(1.0, 1.0 - eigen.eigenValues.x())),
                           0.0,
                           0.0,
                           0.0,
                           std::max(1.0e-4, std::min(1.0, 1.0 - eigen.eigenValues.y())),
                           0.0,
                           0.0,
                           0.0,
                           std::max(1.0e-4, std::min(1.0, 1.0 - eigen.eigenValues.z())));
  result.rotationalTransform(eigen.eigenVectors);
  ENSURE(result.eigenValues().minElement() >= 0.0 &&
         result.eigenValues().maxElement() <= 1.0);
  return result;
}

}

#endif
