#ifndef __PBGWRAPS_SOLIDSPHTYPES__
#define __PBGWRAPS_SOLIDSPHTYPES__

#include "Geometry/Dimension.hh"
#include "SolidSPH/SolidSPHHydroBase.hh"
#include "SolidSPH/NodeCoupling.hh"
#include "SolidSPH/DamagedNodeCoupling.hh"

namespace Spheral {

typedef DamagedNodeCoupling<Dim<1> > DamagedNodeCoupling1d;
typedef DamagedNodeCoupling<Dim<2> > DamagedNodeCoupling2d;
typedef DamagedNodeCoupling<Dim<3> > DamagedNodeCoupling3d;

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
