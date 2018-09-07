#ifndef __PBGWRAPS_HELMHOLTZTYPES__
#define __PBGWRAPS_HELMHOLTZTYPES__

#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/EquationOfState.hh"
#include "Material/HelmholtzEquationOfState.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef HelmholtzEquationOfState<Dim<1> > HelmholtzEquationOfState1d;
typedef HelmholtzEquationOfState<Dim<2> > HelmholtzEquationOfState2d;
typedef HelmholtzEquationOfState<Dim<3> > HelmholtzEquationOfState3d;

}

#endif
