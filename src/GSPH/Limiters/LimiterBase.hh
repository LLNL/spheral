//---------------------------------Spheral++----------------------------------//
// LimiterBase -- base class for scalar slope limiters
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//

#ifndef __Spheral_LimiterBase_hh__
#define __Spheral_LimiterBase_hh__

namespace Spheral {

template<typename Dimension>
class LimiterBase {

public:

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

  LimiterBase(bool TVD,
              bool symmetric);

  ~LimiterBase();

  virtual
  Scalar slopeLimiter(const Scalar) const ;

  virtual
  Scalar fluxLimiter(const Scalar) const = 0;
  

  virtual bool isSymmetric() const;
  virtual bool isTVD() const;

private:
  bool mTVD;
  bool mSymmetric;
};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SlopeLimiter;
}

#endif

