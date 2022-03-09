#include "Utilities/safeInv.hh"
#include "EinfeldtWaveSpeed.hh"
#include <math.h>

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
EinfeldtWaveSpeed<Dimension>::
EinfeldtWaveSpeed(){}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
EinfeldtWaveSpeed<Dimension>::
~EinfeldtWaveSpeed(){}

//========================================================
// wave speed
//========================================================
template<typename Dimension>
void
EinfeldtWaveSpeed<Dimension>::
waveSpeed(const typename Dimension::Scalar rhoi, 
          const typename Dimension::Scalar rhoj, 
          const typename Dimension::Scalar ci,   
          const typename Dimension::Scalar cj, 
          const typename Dimension::Scalar ui,   
          const typename Dimension::Scalar uj,
                typename Dimension::Scalar& Si,  
                typename Dimension::Scalar& Sj) const {

  const auto tiny = std::numeric_limits<Scalar>::epsilon();
  const auto sRhoi = sqrt(rhoi);
  const auto sRhoj = sqrt(rhoj);
  const auto denom = 1.0/std::max(sRhoi+sRhoj,tiny);
  const auto eta = 0.5*sRhoi*sRhoj*denom*denom;
  const auto d2 = (sRhoi * ci*ci + sRhoj * cj*cj)*denom + eta * (ui-uj)*(ui-uj);
  const auto d = sqrt(d2);

  const auto utilde = (sRhoi * ui + sRhoj * uj)*denom;

  Si = rhoi*(std::max( utilde + d, ui+ci ) - ui);
  Sj = rhoj*(std::min( utilde - d, uj-cj ) - uj);
}


}