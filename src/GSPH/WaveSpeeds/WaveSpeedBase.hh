//---------------------------------Spheral++----------------------------------//
// WaveSpeedBase
//----------------------------------------------------------------------------//
#ifndef __Spheral_WaveSpeedBase_hh__
#define __Spheral_WaveSpeedBase_hh__

namespace Spheral {

template<typename Dimension>
class WaveSpeedBase {

public:

  typedef typename Dimension::Scalar Scalar;

  WaveSpeedBase();

  virtual ~WaveSpeedBase();

  virtual
  void waveSpeed(const Scalar rhoi, 
                 const Scalar rhoj, 
                 const Scalar ci,   
                 const Scalar cj, 
                 const Scalar ui,   
                 const Scalar uj,
                       Scalar& Si,        
                       Scalar& Sj) const = 0;

};


}

#endif

