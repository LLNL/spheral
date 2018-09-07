#ifndef __PBGWRAPS_SPHTYPES__
#define __PBGWRAPS_SPHTYPES__

#include "Geometry/Dimension.hh"
#include "SPH/SPHHydroBase.hh"
#include "SPH/PSPHHydroBase.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "SPH/computeSPHOmegaGradhCorrection.hh"
#include "SPH/SPHHydroBaseRZ.hh"
#include "SPH/SPHHydroBaseGSRZ.hh"
#include "SPH/SolidSPHHydroBase.hh"
#include "SPH/SolidSPHHydroBaseRZ.hh"
#include "SPH/NodeCoupling.hh"
#include "SPH/DamagedNodeCoupling.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef DamagedNodeCoupling<Dim<1> > DamagedNodeCoupling1d;
typedef DamagedNodeCoupling<Dim<2> > DamagedNodeCoupling2d;
typedef DamagedNodeCoupling<Dim<3> > DamagedNodeCoupling3d;

typedef SPHHydroBase<Dim<1> > SPHHydroBase1d;
typedef SPHHydroBase<Dim<2> > SPHHydroBase2d;
typedef SPHHydroBase<Dim<3> > SPHHydroBase3d;

typedef PSPHHydroBase<Dim<1> > PSPHHydroBase1d;
typedef PSPHHydroBase<Dim<2> > PSPHHydroBase2d;
typedef PSPHHydroBase<Dim<3> > PSPHHydroBase3d;

typedef SolidSPHHydroBase<Dim<1> > SolidSPHHydroBase1d;
typedef SolidSPHHydroBase<Dim<2> > SolidSPHHydroBase2d;
typedef SolidSPHHydroBase<Dim<3> > SolidSPHHydroBase3d;

}

#endif
