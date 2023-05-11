//---------------------------------Spheral++----------------------------------//
// RectangularFinitePlane -- solid planar boundary for DEM with finite extent
//                           and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/RectangularFinitePlane.hh"

namespace Spheral {

template<typename Dimension>
RectangularFinitePlane<Dimension>::
RectangularFinitePlane(const Vector& point, const Vector& extent, const Tensor& basis):
  SolidBoundary<Dimension>(),
  mPoint(point),
  mBasis(basis),
  mExtent(extent),
  mVelocity(Vector::zero){
}

template<typename Dimension>
RectangularFinitePlane<Dimension>::
~RectangularFinitePlane(){
}


// needs to get fixed
template<typename Dimension>
typename Dimension::Vector
RectangularFinitePlane<Dimension>::
distance(const Vector& position) const { 
  const auto q = mBasis*(position-mPoint);
  const double signX = (q.x() > 0.0 ? 1. : -1.);
  const double signY = (q.y() > 0.0 ? 1. : -1.);
  const double signZ = (q.z() > 0.0 ? 1. : -1.);
  const auto signedExtent = Vector((q.x() > mExtent.x() ? 1 : 0)*signX*mExtent.x(),
                                   (q.y() > mExtent.y() ? 1 : 0)*signY*mExtent.y(),
                                   (q.z() > mExtent.z() ? 1 : 0)*signZ*mExtent.z());
  // PPPPOOOOOOOOP
  //THIS NEEDS TO BE FIXED.
  return mBasis.Transpose()*(q-signedExtent);

}

template<typename Dimension>
typename Dimension::Vector
RectangularFinitePlane<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
RectangularFinitePlane<Dimension>::
update(const double multiplier, const double t, const double dt) {   
}


}