//---------------------------------Spheral++----------------------------------//
// CircularPlaneSolidBoundary -- solid planar boundary for DEM with finite 
//                               extent and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/CircularPlaneSolidBoundary.hh"

#include <string>
using std::string;

namespace Spheral {

template<typename Dimension>
CircularPlaneSolidBoundary<Dimension>::
CircularPlaneSolidBoundary(const Vector& point, const Vector& normal, const Scalar& extent):
  SolidBoundaryBase<Dimension>(),
  mPoint(point),
  mNormal(normal),
  mExtent(extent),
  mVelocity(Vector::zero){
}

template<typename Dimension>
CircularPlaneSolidBoundary<Dimension>::
~CircularPlaneSolidBoundary(){
}


// needs to get fixed
template<typename Dimension>
typename Dimension::Vector
CircularPlaneSolidBoundary<Dimension>::
distance(const Vector& position) const { 
  const auto p = position - mPoint;
  const Vector pn = p.dot(mNormal)*mNormal;
  const Vector pr0 = p-pn;
  const Vector pr = max(pr0.magnitude()-mExtent,0.0) * pr0.unitVector();
  return pn+pr;

}

template<typename Dimension>
typename Dimension::Vector
CircularPlaneSolidBoundary<Dimension>::
localVelocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
CircularPlaneSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {   
  const auto boundaryKey = "CircularPlaneSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto velocityKey = boundaryKey +"_velocity";
  const auto normalKey = boundaryKey +"_normal";
  state.enroll(pointKey,mPoint);
  state.enroll(pointKey,mVelocity);
  state.enroll(pointKey,mNormal);
}

template<typename Dimension>
void
CircularPlaneSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {
  mPoint += multiplier*mVelocity;
}


//------------------------------------------------------------------------------
// Restart
//------------------------------------------------------------------------------
template<typename Dimension>
void
CircularPlaneSolidBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mPoint, pathName + "/point");
  file.write(mNormal, pathName + "/normal");
  file.write(mExtent, pathName + "/extent");
  file.write(mVelocity, pathName + "/velocity");
}


template<typename Dimension>
void
CircularPlaneSolidBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mPoint, pathName + "/point");
  file.read(mNormal, pathName + "/normal");
  file.read(mExtent, pathName + "/extent");
  file.read(mVelocity, pathName + "/velocity");
}

}