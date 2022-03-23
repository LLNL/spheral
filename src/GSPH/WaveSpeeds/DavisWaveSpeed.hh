//---------------------------------Spheral++----------------------------------//
// DavisWaveSpeed 
//   Davis S., (1988) "Simplified Second-Order Godunov-Type Methods" Siam. Sci. 
//   Stat. Comput., 9(3):445-473
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_DavisWaveSpeed_hh__
#define __Spheral_DavisWaveSpeed_hh__

#include "WaveSpeedBase.hh"

namespace Spheral {

template<typename Dimension>
class DavisWaveSpeed : public WaveSpeedBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  DavisWaveSpeed();

  ~DavisWaveSpeed();

  virtual
  void waveSpeed(const Scalar rhoi, const Scalar rhoj, 
                 const Scalar ci,   const Scalar cj, 
                 const Scalar ui,   const Scalar uj,
                       Scalar& Si,        Scalar& Sj) const override;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DavisWaveSpeed;
}

#endif

