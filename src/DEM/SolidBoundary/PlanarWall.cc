//---------------------------------Spheral++----------------------------------//
// SolidBoundary -- this is a base class for solid wall contact bcs
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
  mNormal(normal){
}

template<typename Dimension>
PlanarWall<Dimension>::
~PlanarWall(){
}

template<typename Dimension>
typename Dimension::Scalar
PlanarWall<Dimension>::
value(const Vector& position) const { 
  return (position - mPoint).dot(mNormal);
}

template<typename Dimension>
void
PlanarWall<Dimension>::
update(const double multiplier, const double t, const double dt) {   
}


}