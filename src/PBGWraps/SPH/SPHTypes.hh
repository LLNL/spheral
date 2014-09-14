#ifndef __PBGWRAPS_SPHTYPES__
#define __PBGWRAPS_SPHTYPES__

#include "Geometry/Dimension.hh"
#include "SPH/SPHHydroBase.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "SPH/computeSPHOmegaGradhCorrection.hh"

namespace Spheral {
namespace SPHSpace {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef SPHHydroBase<Dim<1> > SPHHydroBase1d;
typedef SPHHydroBase<Dim<2> > SPHHydroBase2d;
typedef SPHHydroBase<Dim<3> > SPHHydroBase3d;

}
}

#endif
