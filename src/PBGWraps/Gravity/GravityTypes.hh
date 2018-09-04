#ifndef __PBGWRAPS_GRAVITYTYPES__
#define __PBGWRAPS_GRAVITYTYPES__

#include "Geometry/Dimension.hh"
#include "Gravity/NBodyGravity.hh"
#include "Gravity/TreeGravity.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef NBodyGravity<Dim<1> > NBodyGravity1d;
typedef NBodyGravity<Dim<2> > NBodyGravity2d;
typedef NBodyGravity<Dim<3> > NBodyGravity3d;

typedef TreeGravity<Dim<2> > QuadTreeGravity;
typedef TreeGravity<Dim<3> > OctTreeGravity;

}

#endif
