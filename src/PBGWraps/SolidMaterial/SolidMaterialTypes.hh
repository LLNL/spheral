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

typedef LinearPolynomialEquationOfState<Dim<1> > LinearPolynomialEquationOfState1d;
typedef LinearPolynomialEquationOfState<Dim<2> > LinearPolynomialEquationOfState2d;
typedef LinearPolynomialEquationOfState<Dim<3> > LinearPolynomialEquationOfState3d;

typedef GruneisenEquationOfState<Dim<1> > GruneisenEquationOfState1d;
typedef GruneisenEquationOfState<Dim<2> > GruneisenEquationOfState2d;
typedef GruneisenEquationOfState<Dim<3> > GruneisenEquationOfState3d;

typedef TillotsonEquationOfState<Dim<1> > TillotsonEquationOfState1d;
typedef TillotsonEquationOfState<Dim<2> > TillotsonEquationOfState2d;
typedef TillotsonEquationOfState<Dim<3> > TillotsonEquationOfState3d;

typedef MurnahanEquationOfState<Dim<1> > MurnahanEquationOfState1d;
typedef MurnahanEquationOfState<Dim<2> > MurnahanEquationOfState2d;
typedef MurnahanEquationOfState<Dim<3> > MurnahanEquationOfState3d;

typedef SteinbergGuinanStrength<Dim<1> > SteinbergGuinanStrength1d;
typedef SteinbergGuinanStrength<Dim<2> > SteinbergGuinanStrength2d;
typedef SteinbergGuinanStrength<Dim<3> > SteinbergGuinanStrength3d;

typedef SteinbergGuinanLundStrength<Dim<1> > SteinbergGuinanLundStrength1d;
typedef SteinbergGuinanLundStrength<Dim<2> > SteinbergGuinanLundStrength2d;
typedef SteinbergGuinanLundStrength<Dim<3> > SteinbergGuinanLundStrength3d;

}
}

#endif
