#ifndef __PBGWRAPS_HYDROTYPES__
#define __PBGWRAPS_HYDROTYPES__

#include "Geometry/Dimension.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/SecondMomentHourglassControl.hh"
#include "Hydro/ThirdMomentHourglassControl.hh"
#include "Hydro/VoronoiHourglassControl.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef SecondMomentHourglassControl<Dim<1> > SecondMomentHourglassControl1d;
typedef SecondMomentHourglassControl<Dim<2> > SecondMomentHourglassControl2d;
typedef SecondMomentHourglassControl<Dim<3> > SecondMomentHourglassControl3d;

typedef ThirdMomentHourglassControl<Dim<1> > ThirdMomentHourglassControl1d;
typedef ThirdMomentHourglassControl<Dim<2> > ThirdMomentHourglassControl2d;
typedef ThirdMomentHourglassControl<Dim<3> > ThirdMomentHourglassControl3d;

typedef VoronoiHourglassControl<Dim<1> > VoronoiHourglassControl1d;
typedef VoronoiHourglassControl<Dim<2> > VoronoiHourglassControl2d;
typedef VoronoiHourglassControl<Dim<3> > VoronoiHourglassControl3d;

}

#endif
