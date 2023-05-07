//---------------------------------Spheral++----------------------------------//
// InfinitePlane -- rigid planar wall contact for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/InfinitePlane.hh"

namespace Spheral {

template<typename Dimension>
InfinitePlane<Dimension>::
InfinitePlane(const Vector& point, const Vector& normal):
  SolidBoundary<Dimension>(),
  mPoint(point),
  mNormal(normal),
  mVelocity(Vector::zero){
}

template<typename Dimension>
InfinitePlane<Dimension>::
~InfinitePlane(){
}

template<typename Dimension>
typename Dimension::Vector
InfinitePlane<Dimension>::
distance(const Vector& position) const { 
  return (position - mPoint).dot(mNormal)*mNormal;
}

template<typename Dimension>
typename Dimension::Vector
InfinitePlane<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
InfinitePlane<Dimension>::
update(const double multiplier, const double t, const double dt) {   
}


}