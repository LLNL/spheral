//---------------------------------Spheral++----------------------------------//
// RectangularPlaneSolidBoundary -- solid planar boundary for DEM with finite 
//                                  extent and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DEM/SolidBoundary/RectangularPlaneSolidBoundary.hh"

#include <string>
using std::string;

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
  const auto q0 = elementWiseMax(elementWiseMin(q,mExtent),-mExtent);
  return mBasis.Transpose()*(q-q0);
}

template<typename Dimension>
typename Dimension::Vector
RectangularPlaneSolidBoundary<Dimension>::
localVelocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
RectangularPlaneSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {   
  const auto boundaryKey = "RectangularPlaneSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto velocityKey = boundaryKey +"_velocity";
  state.enroll(pointKey,mPoint);
  state.enroll(velocityKey,mVelocity);
}
template<typename Dimension>
void
RectangularPlaneSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
}

//------------------------------------------------------------------------------
// Restart
//------------------------------------------------------------------------------
template<typename Dimension>
void
RectangularPlaneSolidBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mPoint, pathName + "/point");
  file.write(mBasis, pathName + "/basis");
  file.write(mExtent, pathName + "/extent");
  file.write(mVelocity, pathName + "/velocity");
}


template<typename Dimension>
void
RectangularPlaneSolidBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mPoint, pathName + "/point");
  file.read(mBasis, pathName + "/basis");
  file.read(mExtent, pathName + "/extent");
  file.read(mVelocity, pathName + "/velocity");
}


}