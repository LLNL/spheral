//---------------------------------Spheral++----------------------------------//
// VanLeerLimiter
//   Van Leer, B. (1979), "Towards the ultimate conservative difference scheme 
//   V. A second order sequel to Godunov's method", J. Comput. Phys., 32 (1): 
//   101–136
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_VanLeerLimiter_hh__
#define __Spheral_VanLeerLimiter_hh__

#include "LimiterBase.hh"

namespace Spheral {

template<typename Dimension>
class VanLeerLimiter : public LimiterBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  VanLeerLimiter();

  ~VanLeerLimiter();

  virtual
  Scalar fluxLimiter(const Scalar) const override;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class VanLeerLimiter;
}

#endif

