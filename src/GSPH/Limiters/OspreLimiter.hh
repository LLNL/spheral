//---------------------------------Spheral++----------------------------------//
// OspreLimiter 
//   Waterson, N.P.; Deconinck, H. (1995), "A unified approach to the design 
//   and application of bounded higher-order convection schemes," Journal
//   of Computational Physics, 224 (1): 182-207. 
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_OspreLimiter_hh__
#define __Spheral_OspreLimiter_hh__

#include "LimiterBase.hh"

namespace Spheral {

template<typename Dimension>
class OspreLimiter : public LimiterBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  OspreLimiter();

  ~OspreLimiter();

  virtual
  Scalar fluxLimiter(const Scalar) const override;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class OspreLimiter;
}

#endif

