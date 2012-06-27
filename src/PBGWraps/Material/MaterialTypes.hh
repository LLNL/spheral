#ifndef __PBGWRAPS_MATERIALTYPES__
#define __PBGWRAPS_MATERIALTYPES__

#include "Geometry/Dimension.hh"
#include "Material/EquationOfState.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/CGSUnits.hh"
#include "Material/MKSUnits.hh"
#include "Material/CosmologicalUnits.hh"
#include "Material/SolarUnits.hh"
#include "Material/GammaLawGas.hh"
#include "Material/PolytropicEquationOfState.hh"
#include "Material/IsothermalEquationOfState.hh"

namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef PhysicalConstants<CGSUnits> CGS;
typedef PhysicalConstants<MKSUnits> MKS;
typedef PhysicalConstants<CosmologicalUnits> Cosmological;
typedef PhysicalConstants<SolarUnits> Solar;

typedef EquationOfState<Dim<1> > EquationOfState1d;
typedef EquationOfState<Dim<2> > EquationOfState2d;
typedef EquationOfState<Dim<3> > EquationOfState3d;

// CGS
typedef GammaLawGas<Dim<1>, CGS> GammaLawGasCGS1d;
typedef GammaLawGas<Dim<2>, CGS> GammaLawGasCGS2d;
typedef GammaLawGas<Dim<3>, CGS> GammaLawGasCGS3d;

typedef PolytropicEquationOfState<Dim<1>, CGS> PolytropicEquationOfStateCGS1d;
typedef PolytropicEquationOfState<Dim<2>, CGS> PolytropicEquationOfStateCGS2d;
typedef PolytropicEquationOfState<Dim<3>, CGS> PolytropicEquationOfStateCGS3d;

typedef IsothermalEquationOfState<Dim<1>, CGS> IsothermalEquationOfStateCGS1d;
typedef IsothermalEquationOfState<Dim<2>, CGS> IsothermalEquationOfStateCGS2d;
typedef IsothermalEquationOfState<Dim<3>, CGS> IsothermalEquationOfStateCGS3d;

// MKS
typedef GammaLawGas<Dim<1>, MKS> GammaLawGasMKS1d;
typedef GammaLawGas<Dim<2>, MKS> GammaLawGasMKS2d;
typedef GammaLawGas<Dim<3>, MKS> GammaLawGasMKS3d;

typedef PolytropicEquationOfState<Dim<1>, MKS> PolytropicEquationOfStateMKS1d;
typedef PolytropicEquationOfState<Dim<2>, MKS> PolytropicEquationOfStateMKS2d;
typedef PolytropicEquationOfState<Dim<3>, MKS> PolytropicEquationOfStateMKS3d;

typedef IsothermalEquationOfState<Dim<1>, MKS> IsothermalEquationOfStateMKS1d;
typedef IsothermalEquationOfState<Dim<2>, MKS> IsothermalEquationOfStateMKS2d;
typedef IsothermalEquationOfState<Dim<3>, MKS> IsothermalEquationOfStateMKS3d;

// Cosmological
typedef GammaLawGas<Dim<1>, Cosmological> GammaLawGasCosmological1d;
typedef GammaLawGas<Dim<2>, Cosmological> GammaLawGasCosmological2d;
typedef GammaLawGas<Dim<3>, Cosmological> GammaLawGasCosmological3d;

typedef PolytropicEquationOfState<Dim<1>, Cosmological> PolytropicEquationOfStateCosmological1d;
typedef PolytropicEquationOfState<Dim<2>, Cosmological> PolytropicEquationOfStateCosmological2d;
typedef PolytropicEquationOfState<Dim<3>, Cosmological> PolytropicEquationOfStateCosmological3d;

typedef IsothermalEquationOfState<Dim<1>, Cosmological> IsothermalEquationOfStateCosmological1d;
typedef IsothermalEquationOfState<Dim<2>, Cosmological> IsothermalEquationOfStateCosmological2d;
typedef IsothermalEquationOfState<Dim<3>, Cosmological> IsothermalEquationOfStateCosmological3d;

// Solar
typedef GammaLawGas<Dim<1>, Solar> GammaLawGasSolar1d;
typedef GammaLawGas<Dim<2>, Solar> GammaLawGasSolar2d;
typedef GammaLawGas<Dim<3>, Solar> GammaLawGasSolar3d;

typedef PolytropicEquationOfState<Dim<1>, Solar> PolytropicEquationOfStateSolar1d;
typedef PolytropicEquationOfState<Dim<2>, Solar> PolytropicEquationOfStateSolar2d;
typedef PolytropicEquationOfState<Dim<3>, Solar> PolytropicEquationOfStateSolar3d;

typedef IsothermalEquationOfState<Dim<1>, Solar> IsothermalEquationOfStateSolar1d;
typedef IsothermalEquationOfState<Dim<2>, Solar> IsothermalEquationOfStateSolar2d;
typedef IsothermalEquationOfState<Dim<3>, Solar> IsothermalEquationOfStateSolar3d;

}
}

#endif
