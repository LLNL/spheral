//---------------------------------Spheral++----------------------------------//
// AcousticWaveSpeed -- simple acoustic wave speed
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_AcousticWaveSpeed_hh__
#define __Spheral_AcousticWaveSpeed_hh__

#include "WaveSpeedBase.hh"

namespace Spheral {

template<typename Dimension>
class AcousticWaveSpeed : public WaveSpeedBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  AcousticWaveSpeed();

  ~AcousticWaveSpeed();

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
  template<typename Dimension> class AcousticWaveSpeed;
}

#endif

