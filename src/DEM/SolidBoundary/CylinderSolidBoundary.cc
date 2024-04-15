//---------------------------------Spheral++----------------------------------//
// CylinderSolidBoundary -- cylinder with finite length solid boundary for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/CylinderSolidBoundary.hh"

namespace Spheral {

template<typename Dimension>
CylinderSolidBoundary<Dimension>::
CylinderSolidBoundary(const Vector& point, 
               const Vector& axis, 
               const Scalar radius, 
               const Scalar length):
  SolidBoundaryBase<Dimension>(),
  mPoint(point),
  mAxis(axis),
  mRadius(radius),
  mLength(length),
  mVelocity(Vector::zero){
}

template<typename Dimension>
CylinderSolidBoundary<Dimension>::
~CylinderSolidBoundary(){
}

template<typename Dimension>
typename Dimension::Vector
CylinderSolidBoundary<Dimension>::
distance(const Vector& position) const { 
  const auto p = position-mPoint;
  const auto pnMag = p.dot(mAxis);
  const auto pn = pnMag * mAxis;
  const auto paxis = (pnMag > 0 ? max(pnMag-mLength,0.0) : pnMag)*mAxis;
  const auto pr = p-pn;
  return (pr.magnitude() - mRadius)*pr.unitVector() + paxis;
}

template<typename Dimension>
typename Dimension::Vector
CylinderSolidBoundary<Dimension>::
localVelocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {   
  const auto boundaryKey = "CylinderSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto velocityKey = boundaryKey +"_velocity";
  //const auto normalKey = boundaryKey +"_normal";
  state.enrollAny(pointKey,mPoint);
  state.enrollAny(pointKey,mVelocity);
  //state.enrollAny(pointKey,mNormal);
}

template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
}


}