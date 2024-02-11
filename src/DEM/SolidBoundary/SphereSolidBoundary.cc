//---------------------------------Spheral++----------------------------------//
// SphereSolidBoundary -- N-dimensional spherical solid boundary for DEM.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/SphereSolidBoundary.hh"

#include <cmath>

namespace Spheral {

template<typename Dimension>
SphereSolidBoundary<Dimension>::
SphereSolidBoundary(const Vector& center,
                    const Scalar  radius,
                    const Vector& clipPoint,
                    const Vector& clipAxis):
  SolidBoundaryBase<Dimension>(),
  mCenter(center),
  mRadius(radius),
  mClipPoint(clipPoint),
  mClipAxis(clipAxis),
  mClipIntersectionRadius(0.0),
  mVelocity(Vector::zero){
    this->setClipIntersectionRadius();
    mClipAxis = mClipAxis.unitVector();
}

template<typename Dimension>
SphereSolidBoundary<Dimension>::
~SphereSolidBoundary(){
}


template<typename Dimension>
typename Dimension::Vector
SphereSolidBoundary<Dimension>::
distance(const Vector& position) const { 

  // contacting sphere
  const auto contactPoint = (position - mCenter).unitVector()*mRadius + mCenter;
  Vector dist = position - contactPoint;

  const auto planeSignedDistance = (contactPoint - mClipPoint).dot(mClipAxis);

  // if contant pt above clip plane check for edge contact
  if (planeSignedDistance > 0.0){

    // break into perp and in-plane components
    const auto q = (position-mClipPoint);
    const auto qnMag =  q.dot(mClipAxis);
    const auto qn = qnMag * mClipAxis;
    const auto qr = q - qn;

    // if outside circle enforce planar solid bc
    dist = min(qr.magnitude() - mClipIntersectionRadius,0.0)*qr.unitVector() + qn;

  }
  return dist;
}

template<typename Dimension>
typename Dimension::Vector
SphereSolidBoundary<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
SphereSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state,
              const std::string& boundaryKey) {   
}

template<typename Dimension>
void
SphereSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mCenter += multiplier*mVelocity;
}


template<typename Dimension>
void
SphereSolidBoundary<Dimension>::
setClipIntersectionRadius() {
  const auto rcMag = (mClipPoint - mCenter).dot(mClipAxis);
  mClipIntersectionRadius = (rcMag < mRadius ? std::sqrt(mRadius*mRadius-rcMag*rcMag) : 0.0);
  mClipPoint = rcMag * mClipAxis + mCenter;
}

}