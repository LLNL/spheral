#ifndef __PBGWRAPS_SOLIDMATERIALTYPES__
#define __PBGWRAPS_SOLIDMATERIALTYPES__

#include "Geometry/Dimension.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "SolidMaterial/LinearPolynomialEquationOfState.hh"
#include "SolidMaterial/GruneisenEquationOfState.hh"
#include "SolidMaterial/TillotsonEquationOfState.hh"
#include "SolidMaterial/MurnahanEquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "SolidMaterial/ConstantStrength.hh"
#include "SolidMaterial/NullStrength.hh"
#include "SolidMaterial/PolynomialFit.hh"
#include "SolidMaterial/SteinbergGuinanStrength.hh"
#include "SolidMaterial/SteinbergGuinanLundStrength.hh"
#include "SolidMaterial/PorousEquationOfState.hh"
#include "SolidMaterial/StrainPorosity.hh"

using namespace Spheral::Material;

namespace Spheral {
namespace SolidMaterial {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef SolidEquationOfState<Dim<1> > SolidEquationOfState1d;
typedef SolidEquationOfState<Dim<2> > SolidEquationOfState2d;
typedef SolidEquationOfState<Dim<3> > SolidEquationOfState3d;

typedef PorousEquationOfState<Dim<1> > PorousEquationOfState1d;
typedef PorousEquationOfState<Dim<2> > PorousEquationOfState2d;
typedef PorousEquationOfState<Dim<3> > PorousEquationOfState3d;

typedef StrainPorosity<Dim<1> > StrainPorosity1d;
typedef StrainPorosity<Dim<2> > StrainPorosity2d;
typedef StrainPorosity<Dim<3> > StrainPorosity3d;

// CGS
typedef LinearPolynomialEquationOfState<Dim<1>, CGS> LinearPolynomialEquationOfStateCGS1d;
typedef LinearPolynomialEquationOfState<Dim<2>, CGS> LinearPolynomialEquationOfStateCGS2d;
typedef LinearPolynomialEquationOfState<Dim<3>, CGS> LinearPolynomialEquationOfStateCGS3d;

typedef GruneisenEquationOfState<Dim<1>, CGS> GruneisenEquationOfStateCGS1d;
typedef GruneisenEquationOfState<Dim<2>, CGS> GruneisenEquationOfStateCGS2d;
typedef GruneisenEquationOfState<Dim<3>, CGS> GruneisenEquationOfStateCGS3d;

typedef TillotsonEquationOfState<Dim<1>, CGS> TillotsonEquationOfStateCGS1d;
typedef TillotsonEquationOfState<Dim<2>, CGS> TillotsonEquationOfStateCGS2d;
typedef TillotsonEquationOfState<Dim<3>, CGS> TillotsonEquationOfStateCGS3d;

typedef MurnahanEquationOfState<Dim<1>, CGS> MurnahanEquationOfStateCGS1d;
typedef MurnahanEquationOfState<Dim<2>, CGS> MurnahanEquationOfStateCGS2d;
typedef MurnahanEquationOfState<Dim<3>, CGS> MurnahanEquationOfStateCGS3d;

typedef SteinbergGuinanStrength<Dim<1>, CGS> SteinbergGuinanStrengthCGS1d;
typedef SteinbergGuinanStrength<Dim<2>, CGS> SteinbergGuinanStrengthCGS2d;
typedef SteinbergGuinanStrength<Dim<3>, CGS> SteinbergGuinanStrengthCGS3d;

typedef SteinbergGuinanLundStrength<Dim<1>, CGS> SteinbergGuinanLundStrengthCGS1d;
typedef SteinbergGuinanLundStrength<Dim<2>, CGS> SteinbergGuinanLundStrengthCGS2d;
typedef SteinbergGuinanLundStrength<Dim<3>, CGS> SteinbergGuinanLundStrengthCGS3d;

// MKS
typedef LinearPolynomialEquationOfState<Dim<1>, MKS> LinearPolynomialEquationOfStateMKS1d;
typedef LinearPolynomialEquationOfState<Dim<2>, MKS> LinearPolynomialEquationOfStateMKS2d;
typedef LinearPolynomialEquationOfState<Dim<3>, MKS> LinearPolynomialEquationOfStateMKS3d;

typedef GruneisenEquationOfState<Dim<1>, MKS> GruneisenEquationOfStateMKS1d;
typedef GruneisenEquationOfState<Dim<2>, MKS> GruneisenEquationOfStateMKS2d;
typedef GruneisenEquationOfState<Dim<3>, MKS> GruneisenEquationOfStateMKS3d;

typedef TillotsonEquationOfState<Dim<1>, MKS> TillotsonEquationOfStateMKS1d;
typedef TillotsonEquationOfState<Dim<2>, MKS> TillotsonEquationOfStateMKS2d;
typedef TillotsonEquationOfState<Dim<3>, MKS> TillotsonEquationOfStateMKS3d;

typedef MurnahanEquationOfState<Dim<1>, MKS> MurnahanEquationOfStateMKS1d;
typedef MurnahanEquationOfState<Dim<2>, MKS> MurnahanEquationOfStateMKS2d;
typedef MurnahanEquationOfState<Dim<3>, MKS> MurnahanEquationOfStateMKS3d;

typedef SteinbergGuinanStrength<Dim<1>, MKS> SteinbergGuinanStrengthMKS1d;
typedef SteinbergGuinanStrength<Dim<2>, MKS> SteinbergGuinanStrengthMKS2d;
typedef SteinbergGuinanStrength<Dim<3>, MKS> SteinbergGuinanStrengthMKS3d;

typedef SteinbergGuinanLundStrength<Dim<1>, MKS> SteinbergGuinanLundStrengthMKS1d;
typedef SteinbergGuinanLundStrength<Dim<2>, MKS> SteinbergGuinanLundStrengthMKS2d;
typedef SteinbergGuinanLundStrength<Dim<3>, MKS> SteinbergGuinanLundStrengthMKS3d;

// Solar
typedef LinearPolynomialEquationOfState<Dim<1>, Solar> LinearPolynomialEquationOfStateSolar1d;
typedef LinearPolynomialEquationOfState<Dim<2>, Solar> LinearPolynomialEquationOfStateSolar2d;
typedef LinearPolynomialEquationOfState<Dim<3>, Solar> LinearPolynomialEquationOfStateSolar3d;

typedef GruneisenEquationOfState<Dim<1>, Solar> GruneisenEquationOfStateSolar1d;
typedef GruneisenEquationOfState<Dim<2>, Solar> GruneisenEquationOfStateSolar2d;
typedef GruneisenEquationOfState<Dim<3>, Solar> GruneisenEquationOfStateSolar3d;

typedef TillotsonEquationOfState<Dim<1>, Solar> TillotsonEquationOfStateSolar1d;
typedef TillotsonEquationOfState<Dim<2>, Solar> TillotsonEquationOfStateSolar2d;
typedef TillotsonEquationOfState<Dim<3>, Solar> TillotsonEquationOfStateSolar3d;

typedef MurnahanEquationOfState<Dim<1>, Solar> MurnahanEquationOfStateSolar1d;
typedef MurnahanEquationOfState<Dim<2>, Solar> MurnahanEquationOfStateSolar2d;
typedef MurnahanEquationOfState<Dim<3>, Solar> MurnahanEquationOfStateSolar3d;

typedef SteinbergGuinanStrength<Dim<1>, Solar> SteinbergGuinanStrengthSolar1d;
typedef SteinbergGuinanStrength<Dim<2>, Solar> SteinbergGuinanStrengthSolar2d;
typedef SteinbergGuinanStrength<Dim<3>, Solar> SteinbergGuinanStrengthSolar3d;

typedef SteinbergGuinanLundStrength<Dim<1>, Solar> SteinbergGuinanLundStrengthSolar1d;
typedef SteinbergGuinanLundStrength<Dim<2>, Solar> SteinbergGuinanLundStrengthSolar2d;
typedef SteinbergGuinanLundStrength<Dim<3>, Solar> SteinbergGuinanLundStrengthSolar3d;

}
}

#endif
