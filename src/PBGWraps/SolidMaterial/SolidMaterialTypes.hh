#ifndef __PBGWRAPS_SOLIDMATERIALTYPES__
#define __PBGWRAPS_SOLIDMATERIALTYPES__

#include "Geometry/Dimension.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "SolidMaterial/LinearPolynomialEquationOfState.hh"
#include "SolidMaterial/GruneisenEquationOfState.hh"
#include "SolidMaterial/OsborneEquationOfState.hh"
#include "SolidMaterial/TillotsonEquationOfState.hh"
#include "SolidMaterial/MurnahanEquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "SolidMaterial/ConstantStrength.hh"
#include "SolidMaterial/NullStrength.hh"
#include "SolidMaterial/PolynomialFit.hh"
#include "SolidMaterial/SteinbergGuinanStrength.hh"
#include "SolidMaterial/SteinbergGuinanLundStrength.hh"
#include "SolidMaterial/JohnsonCookStrength.hh"
#include "SolidMaterial/CollinsStrength.hh"
#include "SolidMaterial/PorousEquationOfState.hh"
#include "SolidMaterial/PorousStrengthModel.hh"
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

typedef OsborneEquationOfState<Dim<1> > OsborneEquationOfState1d;
typedef OsborneEquationOfState<Dim<2> > OsborneEquationOfState2d;
typedef OsborneEquationOfState<Dim<3> > OsborneEquationOfState3d;

typedef TillotsonEquationOfState<Dim<1> > TillotsonEquationOfState1d;
typedef TillotsonEquationOfState<Dim<2> > TillotsonEquationOfState2d;
typedef TillotsonEquationOfState<Dim<3> > TillotsonEquationOfState3d;

typedef MurnahanEquationOfState<Dim<1> > MurnahanEquationOfState1d;
typedef MurnahanEquationOfState<Dim<2> > MurnahanEquationOfState2d;
typedef MurnahanEquationOfState<Dim<3> > MurnahanEquationOfState3d;

typedef StrengthModel<Dim<1> > StrengthModel1d;
typedef StrengthModel<Dim<2> > StrengthModel2d;
typedef StrengthModel<Dim<3> > StrengthModel3d;

typedef NullStrength<Dim<1> > NullStrength1d;
typedef NullStrength<Dim<2> > NullStrength2d;
typedef NullStrength<Dim<3> > NullStrength3d;

typedef ConstantStrength<Dim<1> > ConstantStrength1d;
typedef ConstantStrength<Dim<2> > ConstantStrength2d;
typedef ConstantStrength<Dim<3> > ConstantStrength3d;

typedef SteinbergGuinanStrength<Dim<1> > SteinbergGuinanStrength1d;
typedef SteinbergGuinanStrength<Dim<2> > SteinbergGuinanStrength2d;
typedef SteinbergGuinanStrength<Dim<3> > SteinbergGuinanStrength3d;

typedef SteinbergGuinanLundStrength<Dim<1> > SteinbergGuinanLundStrength1d;
typedef SteinbergGuinanLundStrength<Dim<2> > SteinbergGuinanLundStrength2d;
typedef SteinbergGuinanLundStrength<Dim<3> > SteinbergGuinanLundStrength3d;

typedef JohnsonCookStrength<Dim<1> > JohnsonCookStrength1d;
typedef JohnsonCookStrength<Dim<2> > JohnsonCookStrength2d;
typedef JohnsonCookStrength<Dim<3> > JohnsonCookStrength3d;

typedef CollinsStrength<Dim<1> > CollinsStrength1d;
typedef CollinsStrength<Dim<2> > CollinsStrength2d;
typedef CollinsStrength<Dim<3> > CollinsStrength3d;

typedef PorousStrengthModel<Dim<1> > PorousStrengthModel1d;
typedef PorousStrengthModel<Dim<2> > PorousStrengthModel2d;
typedef PorousStrengthModel<Dim<3> > PorousStrengthModel3d;

}
}

#endif
