#ifndef __PBGWRAPS_GRAVITYTYPES__
#define __PBGWRAPS_GRAVITYTYPES__

#include "Geometry/Dimension.hh"
#include "Gravity/NBodyGravity.hh"

namespace Spheral {
namespace GravitySpace {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef NBodyGravity<Dim<1> > NBodyGravity1d;
typedef NBodyGravity<Dim<2> > NBodyGravity2d;
typedef NBodyGravity<Dim<3> > NBodyGravity3d;

}
}

#endif
