//------------------------------------------------------------------------------
// The entropy weighted energy form.
// Useful for apportioning work between points -- used in 
// SpecificThermalEnergyPolicy for instance.
//------------------------------------------------------------------------------
#ifndef __Spheral_entropyWeightingFunction__
#define __Spheral_entropyWeightingFunction__

#include "Utilities/DBC.hh"

#include <algorithm>

namespace Spheral {

inline
double
entropyWeighting(const double si,
                 const double sj,
                 const double duij) {
  double result = 0.5;
  const double smin = std::min(abs(si), abs(sj));
  const double smax = std::max(abs(si), abs(sj));
  if (smax > 1.0e-15) {
    CHECK(smin + smax > 1.0e-15);
    if (duij > 0.0) {    // Heating
      if (si > sj) {
        result = smin/(smin + smax);
      } else {
        result = smax/(smin + smax);
      }
    } else {             // Cooling
      if (si > sj) {
        result = smax/(smin + smax);
      } else {
        result = smin/(smin + smax);
      }
    }
  }
  CHECK(result >= 0.0 and result <= 1.0);
  return result;
}

}

#endif
