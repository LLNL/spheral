//---------------------------------Spheral++----------------------------------//
// CircularFinitePlane -- solid planar boundary for DEM with finite extent
//                           and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/CircularFinitePlane.hh"

namespace Spheral {

template<typename Dimension>
CircularFinitePlane<Dimension>::
CircularFinitePlane(const Vector& point, const Vector& normal, const Scalar& extent):
  SolidBoundary<Dimension>(),
  mPoint(point),
  mNormal(normal),
  mExtent(extent),
  mVelocity(Vector::zero){
}

template<typename Dimension>
CircularFinitePlane<Dimension>::
~CircularFinitePlane(){
}


// needs to get fixed
template<typename Dimension>
typename Dimension::Vector
CircularFinitePlane<Dimension>::
distance(const Vector& position) const { 
  const auto p = position - mPoint;
  const Vector pn = p.dot(mNormal)*mNormal;
  const Vector pr0 = p-pn;
  const Vector pr = max(pr0.magnitude()-mExtent,0.0) * pr0.unitVector();
  return pn+pr;

}

template<typename Dimension>
typename Dimension::Vector
CircularFinitePlane<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
CircularFinitePlane<Dimension>::
update(const double multiplier, const double t, const double dt) {
  mPoint += multiplier*mVelocity;
}


}