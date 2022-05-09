//---------------------------------Spheral++----------------------------------//
// SuperbeeLimiter 
//   Roe, P.L. (1986), "Characteristic-based schemes for the Euler equations", 
//   Annu. Rev. Fluid Mech., 18: 337–365
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_SuperbeeLimiter_hh__
#define __Spheral_SuperbeeLimiter_hh__

#include "LimiterBase.hh"

namespace Spheral {

template<typename Dimension>
class SuperbeeLimiter : public LimiterBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  SuperbeeLimiter();

  ~SuperbeeLimiter();

  virtual
  Scalar fluxLimiter(const Scalar) const override;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SuperbeeLimiter;
}

#endif

