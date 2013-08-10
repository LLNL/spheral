#ifndef __PBGWRAPS_SVPHTYPES__
#define __PBGWRAPS_SVPHTYPES__

#include "Geometry/Dimension.hh"
#include "SVPH/sampleFieldListSVPH.hh"
#include "SVPH/gradientFieldListSVPH.hh"
#include "SVPH/SVPHHydroBase.hh"

namespace Spheral {
namespace SVPHSpace {

typedef SVPHHydroBase<Dim<1> > SVPHHydroBase1d;
typedef SVPHHydroBase<Dim<2> > SVPHHydroBase2d;
typedef SVPHHydroBase<Dim<3> > SVPHHydroBase3d;

}
}

#endif
