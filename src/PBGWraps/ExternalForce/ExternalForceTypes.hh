#ifndef __PBGWRAPS_EXTERNALFORCETYPES__
#define __PBGWRAPS_EXTERNALFORCETYPES__

#include "Geometry/Dimension.hh"
#include "ExternalForce/PointPotential.hh"
#include "ExternalForce/ConstantAcceleration.hh"
#include "ExternalForce/LinearAcceleration.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef PointPotential<Dim<1> > PointPotential1d;
typedef PointPotential<Dim<2> > PointPotential2d;
typedef PointPotential<Dim<3> > PointPotential3d;

typedef ConstantAcceleration<Dim<1> > ConstantAcceleration1d;
typedef ConstantAcceleration<Dim<2> > ConstantAcceleration2d;
typedef ConstantAcceleration<Dim<3> > ConstantAcceleration3d;

typedef LinearAcceleration<Dim<1> > LinearAcceleration1d;
typedef LinearAcceleration<Dim<2> > LinearAcceleration2d;
typedef LinearAcceleration<Dim<3> > LinearAcceleration3d;

}

#endif
