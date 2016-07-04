#ifndef __Spheral_volumeSpacing__
#define __Spheral_volumeSpacing__

#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the spacing implied by a volume.
//------------------------------------------------------------------------------
template<typename Dimension> double volumeSpacing(const double vol);

template<> 
inline 
double 
volumeSpacing<Dim<1> >(const double vol) {
  return 0.5*vol;
}
  
template<> 
inline
double 
volumeSpacing<Dim<2> >(const double vol) {
  return sqrt(vol/M_PI);
}
  
template<> 
inline
double 
volumeSpacing<Dim<3> >(const double vol) {
  return pow(vol/(4.0/3.0*M_PI), 1.0/3.0);
}

}

#endif
