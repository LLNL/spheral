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
  const double smin = std::min(si, sj);
  const double ssi = si - smin;
  const double ssj = sj - smin;
  CHECK(ssi >= 0.0 and ssj >= 0.0);

  // Work in the shifted positive entropies.  
  double result = 0.5;
  const double ssmin = std::min(ssi, ssj);
  const double ssmax = std::max(ssi, ssj);
  const double sssum = ssmin + ssmax;
  if (sssum > 1.0e-10 and
      ssmax - ssmin > 1.0e-10*sssum) {   // If the entropies are equal, equipartion.

    if (duij > 0.0) {                   // Heating
      if (ssi > ssj) {
        result = ssmin/sssum;
      } else {
        result = ssmax/sssum;
      }

    } else {                            // Cooling
      if (ssi > ssj) {
        result = ssmax/sssum;
      } else {
        result = ssmin/sssum;
      }
    }
  }

  CHECK(result >= 0.0 and result <= 1.0);
  return result;
}

}

#endif
