//---------------------------------Spheral++----------------------------------//
// SlopeLimiterBase -- base class for scalar slope limiters
//----------------------------------------------------------------------------//
#ifndef __Spheral_SlopeLimiterBase_hh__
#define __Spheral_SlopeLimiterBase_hh__

namespace Spheral {

template<typename Dimension>
class SlopeLimiterBase {

public:

  typedef typename Dimension::Scalar Scalar;

  SlopeLimiterBase();

  ~SlopeLimiterBase();

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

