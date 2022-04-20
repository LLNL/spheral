//---------------------------------Spheral++----------------------------------//
// VanAlbaLimiter 
//   Van Albada, G.D.; Van Leer, B.; Roberts, W.W. (1982), 
//   "A comparative study of computational methods in cosmic gas dynamics", 
//   Astronomy and Astrophysics, 108 (1): 76–84
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_VanAlbaLimiter_hh__
#define __Spheral_VanAlbaLimiter_hh__

#include "LimiterBase.hh"

namespace Spheral {

template<typename Dimension>
class VanAlbaLimiter : public LimiterBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  VanAlbaLimiter();

  ~VanAlbaLimiter();

  virtual
  Scalar fluxLimiter(const Scalar) const override;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class VanAlbaLimiter;
}

#endif

