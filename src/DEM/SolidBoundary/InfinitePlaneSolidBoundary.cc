//---------------------------------Spheral++----------------------------------//
// InfinitePlaneSolidBoundary -- rigid planar wall contact for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/InfinitePlaneSolidBoundary.hh"

namespace Spheral {

template<typename Dimension>
InfinitePlaneSolidBoundary<Dimension>::
InfinitePlaneSolidBoundary(const Vector& point, const Vector& normal):
  SolidBoundaryBase<Dimension>(),
  mPoint(point),
  mNormal(normal),
  mVelocity(Vector::zero){
}

template<typename Dimension>
InfinitePlaneSolidBoundary<Dimension>::
~InfinitePlaneSolidBoundary(){
}

template<typename Dimension>
typename Dimension::Vector
InfinitePlaneSolidBoundary<Dimension>::
distance(const Vector& position) const { 
  return (position - mPoint).dot(mNormal)*mNormal;
}

template<typename Dimension>
typename Dimension::Vector
InfinitePlaneSolidBoundary<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
InfinitePlaneSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
}


}