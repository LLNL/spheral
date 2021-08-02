#include "Utilities/safeInv.hh"
#include "EinfeldtWaveSpeed.hh"

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

  const auto sRhoi = FastMath::Sqrt(rhoi);
  const auto sRhoj = FastMath::Sqrt(rhoj);
  const auto denom = safeInv(sRhoi+sRhoj);
  const auto eta = 0.5*sRhoi*sRhoj*denom*denom;
  const auto d2 = (sRhoi * ci*ci + sRhoj * cj*cj)*denom + eta * (ui-uj)*(ui-uj);
  const auto d = FastMath::Sqrt(d2);

  const auto utilde = (sRhoi * ui + sRhoj * uj)*denom;

  Si = rhoi*((utilde + d)-ui);
  Sj = rhoj*((utilde - d)-uj);
}


}