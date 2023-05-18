//---------------------------------Spheral++----------------------------------//
// FiniteCylinder -- solid planar boundary for DEM with finite extent
//                           and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/FiniteCylinder.hh"

namespace Spheral {

template<typename Dimension>
FiniteCylinder<Dimension>::
FiniteCylinder(const Vector& point, 
               const Vector& axis, 
               const Scalar radius, 
               const Scalar length):
  SolidBoundary<Dimension>(),
  mPoint(point),
  mAxis(axis),
  mRadius(radius),
  mLength(length),
  mVelocity(Vector::zero){
}

template<typename Dimension>
FiniteCylinder<Dimension>::
~FiniteCylinder(){
}


template<typename Dimension>
typename Dimension::Vector
FiniteCylinder<Dimension>::
distance(const Vector& position) const { 
  const auto p = position-mPoint;
  const auto pn = p.dot(mAxis);
  const auto pr = p-pn;
  return (pr.magnitude() - mRadius)*pr.unitVector();
}

template<typename Dimension>
typename Dimension::Vector
FiniteCylinder<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
FiniteCylinder<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
}


}