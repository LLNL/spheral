//---------------------------------Spheral++----------------------------------//
// PlanarWall -- rigid planar wall contact for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/PlanarWall.hh"

namespace Spheral {

template<typename Dimension>
PlanarWall<Dimension>::
PlanarWall(const Vector& point, const Vector& normal):
  SolidBoundary<Dimension>(),
  mPoint(point),
  mNormal(normal),
  mVelocity(Vector::zero){
}

template<typename Dimension>
PlanarWall<Dimension>::
~PlanarWall(){
}

template<typename Dimension>
typename Dimension::Vector
PlanarWall<Dimension>::
distance(const Vector& position) const { 
  return (position - mPoint).dot(mNormal)*mNormal;
}

template<typename Dimension>
typename Dimension::Vector
PlanarWall<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
PlanarWall<Dimension>::
update(const double multiplier, const double t, const double dt) {   
}


}