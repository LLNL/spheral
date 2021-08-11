//---------------------------------Spheral++----------------------------------//
// LimiterBase -- base class for scalar slope limiters
//----------------------------------------------------------------------------//
#ifndef __Spheral_LimiterBase_hh__
#define __Spheral_LimiterBase_hh__

namespace Spheral {

template<typename Dimension>
class LimiterBase {

public:

  typedef typename Dimension::Scalar Scalar;

  LimiterBase();

  ~LimiterBase();

  virtual
  Scalar slopeLimiter(const Scalar) const ;

  virtual
  Scalar fluxLimiter(const Scalar) const = 0;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SlopeLimiter;
}

#endif

