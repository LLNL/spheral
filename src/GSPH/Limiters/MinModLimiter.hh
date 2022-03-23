//---------------------------------Spheral++----------------------------------//
// MinModLimiter 
//   Roe, P.L. (1986), "Characteristic-based schemes for the Euler equations", 
//   Annu. Rev. Fluid Mech., 18: 337–365
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_MinModLimiter_hh__
#define __Spheral_MinModLimiter_hh__

#include "LimiterBase.hh"

namespace Spheral {

template<typename Dimension>
class MinModLimiter : public LimiterBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  MinModLimiter();

  ~MinModLimiter();

  virtual
  Scalar fluxLimiter(const Scalar) const override;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class MinModLimiter;
}

#endif

