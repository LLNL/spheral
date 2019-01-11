#ifndef __PBGWRAPS_SVPHTYPES__
#define __PBGWRAPS_SVPHTYPES__

#include "Geometry/Dimension.hh"
#include "SVPH/sampleFieldListSVPH.hh"
#include "SVPH/gradientFieldListSVPH.hh"
#include "SVPH/SVPHHydroBase.hh"
#include "SVPH/SVPHFacetedHydroBase.hh"

namespace Spheral {

typedef SVPHHydroBase<Dim<1> > SVPHHydroBase1d;
typedef SVPHHydroBase<Dim<2> > SVPHHydroBase2d;
typedef SVPHHydroBase<Dim<3> > SVPHHydroBase3d;

typedef SVPHFacetedHydroBase<Dim<1> > SVPHFacetedHydroBase1d;
typedef SVPHFacetedHydroBase<Dim<2> > SVPHFacetedHydroBase2d;
typedef SVPHFacetedHydroBase<Dim<3> > SVPHFacetedHydroBase3d;

}

#endif
