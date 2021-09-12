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
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

  LimiterBase(bool TVD,
              bool symmetric);

  ~LimiterBase();

  virtual
  Scalar slopeLimiter(const Scalar) const ;

  virtual
  Scalar fluxLimiter(const Scalar) const = 0;

  virtual 
  void construct(const Vector& ri,
                 const Vector& rj,
                 const Scalar& yi,
                 const Scalar& yj,
                 const Vector& DyDxi,
                 const Vector& DyDxj,
                       Scalar& ytildei,
                       Scalar& ytildej) const;

  virtual 
  void construct(const Vector& ri,
                 const Vector& rj,
                 const Vector& yi,
                 const Vector& yj,
                 const Tensor& DyDxi,
                 const Tensor& DyDxj,
                       Vector& ytildei,
                       Vector& ytildej) const;
  

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

