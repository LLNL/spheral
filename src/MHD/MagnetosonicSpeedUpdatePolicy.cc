#include "MHD/MagnetosonicSpeedUpdatePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/FluidNodeList.hh"
#include "MHD/MHDFieldNames.hh"
#include "Material/EquationOfState.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
MagnetosonicSpeedUpdatePolicy::
MagnetosonicSpeedUpdatePolicy(double mu0):
   UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >(HydroFieldNames::massDensity,
                                                            HydroFieldNames::specificThermalEnergy,
                                                            MHDFieldNames::magneticInduction),
   mMu0(mu0)
{
} // end constructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
MagnetosonicSpeedUpdatePolicy::
~MagnetosonicSpeedUpdatePolicy()
{
} // end destructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
MagnetosonicSpeedUpdatePolicy::
update(const KeyType& key,
       State<Dim<3> >& state,
       StateDerivatives<Dim<3> >& derivs,
       const double multiplier,
       const double t,
       const double dt)
{
   // The speed of the fast magnetosonic wave is given by
   //
   //  2      2      2
   // C  =  Cs  +  Va
   //
   // where Cs is the sound speed, Va is the Alfven speed, and
   // C is the speed of the fast magnetosonic wave.  The way that we 
   // implement this change is to simply change the sound speed itself:
   //
   // Cs += (C - Cs), where C is the fast magnetosonic speed.
   // Therefore, we'd better be manipulating Cs itself.
   REQUIRE(key.second == HydroFieldNames::soundSpeed);

   // Get the mass density and specific thermal energy fields from the state.
   const KeyType massDensityKey(key.first, HydroFieldNames::massDensity);
   const KeyType energyKey(key.first, HydroFieldNames::specificThermalEnergy);
   CHECK(state.scalarFieldRegistered(massDensityKey));
   CHECK(state.scalarFieldRegistered(energyKey));
   Field<Dim<3>, Dim<3>::Scalar>& soundSpeed = state.scalarField(key);
   const Field<Dim<3>, Dim<3>::Scalar>& massDensity = state.scalarField(massDensityKey);
   const Field<Dim<3>, Dim<3>::Scalar>& energy = state.scalarField(energyKey);

   // Get the eos.  This cast is ugly, but is a work-around for now.
   const FluidNodeList<Dim<3> >* fluidNodeListPtr = 
      (const FluidNodeList<Dim<3> >*) key.first;
   CHECK(fluidNodeListPtr != 0);
   const EquationOfState<Dim<3> >& eos = fluidNodeListPtr->equationOfState();

   // Now set the sound speed.
   eos.setSoundSpeed(soundSpeed, massDensity, energy);

   // Get the magnetic field strength so that we can compute the Alfven 
   // speed.
   const KeyType BKey(key.first, MHDFieldNames::magneticInduction);
   const KeyType rhoKey(key.first, HydroFieldNames::massDensity);
   const Field<Dim<3>, Dim<3>::Vector>& B = state.vectorField(BKey);
   const Field<Dim<3>, Dim<3>::Scalar>& rho = state.scalarField(rhoKey);

   // Go over all the internal nodes and update the sound speed to match the 
   // fast magnetosonic speed.
   for (int i = 0; i < soundSpeed.nodeList().numInternalNodes(); ++i)
   {
      // Compute the square of the fast magnetosonic speed.
      double Va = B[i].magnitude()/sqrt(mMu0*rho[i]);
      double C2 = soundSpeed[i] * soundSpeed[i] + Va * Va;

      // Add the difference to the sound speed.
      double dCs = sqrt(C2) - soundSpeed[i];
      soundSpeed[i] += dCs;
   }

   // We're done.
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
bool
MagnetosonicSpeedUpdatePolicy::
operator==(const Spheral::UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >& rhs) const
{
   // We're only equal if the other guy is the same type.
   const MagnetosonicSpeedUpdatePolicy* rhsPtr = 
      dynamic_cast<const MagnetosonicSpeedUpdatePolicy*>(&rhs);
   if (rhsPtr == 0) {
      return false;
   } else {
      return true;
   }
}
//----------------------------------------------------------------------------

}
