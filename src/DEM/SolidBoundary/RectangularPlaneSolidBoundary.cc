//---------------------------------Spheral++----------------------------------//
// RectangularPlaneSolidBoundary -- solid planar boundary for DEM with finite 
//                                  extent and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/RectangularPlaneSolidBoundary.hh"

namespace Spheral {

template<typename Dimension>
RectangularPlaneSolidBoundary<Dimension>::
RectangularPlaneSolidBoundary(const Vector& point, const Vector& extent, const Tensor& basis):
  SolidBoundaryBase<Dimension>(),
  mPoint(point),
  mBasis(basis),
  mExtent(extent),
  mVelocity(Vector::zero){
}

template<typename Dimension>
RectangularPlaneSolidBoundary<Dimension>::
~RectangularPlaneSolidBoundary(){
}


template<typename Dimension>
typename Dimension::Vector
RectangularPlaneSolidBoundary<Dimension>::
distance(const Vector& position) const { 
  const auto q = mBasis*(position-mPoint);
  const double signX = (q.x() > 0.0 ? 1. : -1.);
  const double signY = (q.y() > 0.0 ? 1. : -1.);
  const double signZ = (q.z() > 0.0 ? 1. : -1.);
  const auto signedExtent = Vector((q.x() > mExtent.x() ? 1 : 0)*signX*mExtent.x(),
                                   (q.y() > mExtent.y() ? 1 : 0)*signY*mExtent.y(),
                                   (q.z() > mExtent.z() ? 1 : 0)*signZ*mExtent.z());
  return mBasis.Transpose()*(q-signedExtent);
}

template<typename Dimension>
typename Dimension::Vector
RectangularPlaneSolidBoundary<Dimension>::
velocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
RectangularPlaneSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
}


}