#ifndef __PBGWRAPS_MATERIALTYPES__
#define __PBGWRAPS_MATERIALTYPES__

#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/EquationOfState.hh"
#include "Material/GammaLawGas.hh"
#include "Material/PolytropicEquationOfState.hh"
#include "Material/IsothermalEquationOfState.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef EquationOfState<Dim<1> > EquationOfState1d;
typedef EquationOfState<Dim<2> > EquationOfState2d;
typedef EquationOfState<Dim<3> > EquationOfState3d;

typedef GammaLawGas<Dim<1> > GammaLawGas1d;
typedef GammaLawGas<Dim<2> > GammaLawGas2d;
typedef GammaLawGas<Dim<3> > GammaLawGas3d;

typedef PolytropicEquationOfState<Dim<1> > PolytropicEquationOfState1d;
typedef PolytropicEquationOfState<Dim<2> > PolytropicEquationOfState2d;
typedef PolytropicEquationOfState<Dim<3> > PolytropicEquationOfState3d;

typedef IsothermalEquationOfState<Dim<1> > IsothermalEquationOfState1d;
typedef IsothermalEquationOfState<Dim<2> > IsothermalEquationOfState2d;
typedef IsothermalEquationOfState<Dim<3> > IsothermalEquationOfState3d;

}

#endif
