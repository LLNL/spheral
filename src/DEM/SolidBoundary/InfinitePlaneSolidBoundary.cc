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
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  const auto boundaryKey = "InfinitePlaneSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto velocityKey = boundaryKey +"_velocity";
  const auto normalKey = boundaryKey +"_normal";
  state.enrollAny(pointKey,mPoint);
  state.enrollAny(velocityKey,mVelocity);
  state.enrollAny(normalKey,mNormal);
}

template<typename Dimension>
void
InfinitePlaneSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
}


}