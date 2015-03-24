#ifndef __PBGWRAPS_SOLIDSPHTYPES__
#define __PBGWRAPS_SOLIDSPHTYPES__

#include "Geometry/Dimension.hh"
#include "SolidSPH/SolidSPHHydroBase.hh"

namespace Spheral {
namespace SolidSPHSpace {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef SolidSPHHydroBase<Dim<1> > SolidSPHHydroBase1d;
typedef SolidSPHHydroBase<Dim<2> > SolidSPHHydroBase2d;
typedef SolidSPHHydroBase<Dim<3> > SolidSPHHydroBase3d;

}
}

#endif
