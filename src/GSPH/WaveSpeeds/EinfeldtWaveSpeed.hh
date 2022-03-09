//---------------------------------Spheral++----------------------------------//
// EinfeldtWaveSpeed
//----------------------------------------------------------------------------//
#ifndef __Spheral_EinfeldtWaveSpeed_hh__
#define __Spheral_EinfeldtWaveSpeed_hh__

#include "WaveSpeedBase.hh"

namespace Spheral {

template<typename Dimension>
class EinfeldtWaveSpeed : public WaveSpeedBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  EinfeldtWaveSpeed();

  ~EinfeldtWaveSpeed();

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
  template<typename Dimension> class EinfeldtWaveSpeed;
}

#endif

