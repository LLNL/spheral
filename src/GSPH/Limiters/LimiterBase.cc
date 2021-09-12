#include "LimiterBase.hh"

namespace Spheral {

//========================================================
// Constructor
//========================================================
template<typename Dimension>
LimiterBase<Dimension>::
LimiterBase(bool TVD,
            bool symmetric):
  mTVD(TVD),
  mSymmetric(symmetric){

}

//========================================================
// Destructor
//========================================================
template<typename Dimension>
LimiterBase<Dimension>::
~LimiterBase(){}

//========================================================
// slope limiter
//========================================================
template<typename Dimension>
typename Dimension::Scalar
LimiterBase<Dimension>::
slopeLimiter(const typename Dimension::Scalar x) const {
  return ( x > 0.0 ? 2.0/(1.0+x) * this->fluxLimiter(x) : 0.0 );
}

//========================================================
// reconstruct from limited gradient
//========================================================
template<typename Dimension>
void
LimiterBase<Dimension>::
construct(const typename Dimension::Vector& ri,
          const typename Dimension::Vector& rj,
          const typename Dimension::Scalar& yi,
          const typename Dimension::Scalar& yj,
          const typename Dimension::Vector& DyDxi,
          const typename Dimension::Vector& DyDxj,
                typename Dimension::Scalar& ytildei,
                typename Dimension::Scalar& ytildej) const {
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto rij = (ri-rj);

  // relavant deltas in field value
  const auto Dy0 = (yi-yj);
  const auto Dyi = 0.5*DyDxi.dot(rij);
  const auto Dyj = 0.5*DyDxj.dot(rij);

  // ratios of SPH derivs to ij particle difference
  const auto denom = 2.0 / (sgn(Dy0) * std::max(tiny,abs(Dy0)));
  const auto ratio0i = Dyi * denom;
  const auto ratio0j = Dyj * denom;

  // limiter function 
  const auto phi0i = this->slopeLimiter(ratio0i);
  const auto phi0j = this->slopeLimiter(ratio0j);
  
  // if our total delta is too large we need reduce it
  const auto ratio0ij = std::min(1.0, 2.0/std::max(tiny,(phi0i*ratio0i+phi0j*ratio0j)));
  const auto phi0ij = this->slopeLimiter(ratio0ij) * ratio0ij;

  // linear constructed inteface values
  ytildei = yi + phi0ij * phi0i * Dyi;
  ytildej = yj + phi0ij * phi0j * Dyj;
}

//========================================================
// reconstruct from limited gradient
//========================================================
template<typename Dimension>
void
LimiterBase<Dimension>::
construct(const typename Dimension::Vector& ri,
          const typename Dimension::Vector& rj,
          const typename Dimension::Vector& yi,
          const typename Dimension::Vector& yj,
          const typename Dimension::Tensor& DyDxi,
          const typename Dimension::Tensor& DyDxj,
                typename Dimension::Vector& ytildei,
                typename Dimension::Vector& ytildej) const {
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto rij = (ri-rj);
  const auto rhatij = rij.unitVector();

  // relavant deltas in field value
  const auto Dy0 = (yi-yj);
  const auto Dyi = 0.5*DyDxi.dot(rij);   
  const auto Dyj = 0.5*DyDxj.dot(rij);  

  // scalar delta along line of action s
  const auto Dy0s = Dy0.dot(rhatij);
  const auto Dyis = Dyi.dot(rhatij);
  const auto Dyjs = Dyj.dot(rhatij);

  // ratios of SPH derivs to ij particle difference
  const auto denom = 2.0 / (sgn(Dy0s) * std::max(tiny,abs(Dy0s)));
  const auto ratio0i = Dyis * denom;
  const auto ratio0j = Dyjs * denom;

  // limiter function 
  const auto phi0i = this->slopeLimiter(ratio0i);
  const auto phi0j = this->slopeLimiter(ratio0j);
  
  const auto ratio0ij = std::min(1.0, 2.0/std::max(tiny,(phi0i*ratio0i+phi0j*ratio0j)));
  const auto phi0ij = this->slopeLimiter(ratio0ij) * ratio0ij;

  ytildei = yi + phi0ij*phi0i * Dyi;
  ytildej = yj + phi0ij*phi0j * Dyj;
}


//========================================================
// is this TVD?
//========================================================
template<typename Dimension>
bool
LimiterBase<Dimension>::
isTVD() const {
  return mTVD;
}

//========================================================
// is this Symmetric?
//========================================================

template<typename Dimension>
bool
LimiterBase<Dimension>::
isSymmetric() const {
  return mSymmetric;
}

}