#ifndef __PBGWRAPS_TaylorSPHTYPES__
#define __PBGWRAPS_TaylorSPHTYPES__

#include "Geometry/Dimension.hh"
#include "TaylorSPH/TaylorSPHHydroBase.hh"
#include "TaylorSPH/computeTaylorSPHCorrections.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef TaylorSPHHydroBase<Dim<1> > TaylorSPHHydroBase1d;
typedef TaylorSPHHydroBase<Dim<2> > TaylorSPHHydroBase2d;
typedef TaylorSPHHydroBase<Dim<3> > TaylorSPHHydroBase3d;

}

#endif
