#include "MHD/SpitzerResistivityUpdatePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "MHD/MHDFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
SpitzerResistivityUpdatePolicy::
SpitzerResistivityUpdatePolicy(double C, double Rmax):
   UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >(HydroFieldNames::massDensity,
                                                            HydroFieldNames::temperature),
   mC(C),
   mRmax(Rmax),
   mRhoMin(1.0e-10)
{
} // end constructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
SpitzerResistivityUpdatePolicy::
SpitzerResistivityUpdatePolicy(double C, double Rmax, double rhoMin):
   UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >(HydroFieldNames::massDensity,
                                                            HydroFieldNames::temperature),
   mC(C),
   mRmax(Rmax),
   mRhoMin(rhoMin)
{
} // end constructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
SpitzerResistivityUpdatePolicy::
~SpitzerResistivityUpdatePolicy()
{
} // end destructor
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
SpitzerResistivityUpdatePolicy::
update(const KeyType& key,
       State<Dim<3> >& state,
       StateDerivatives<Dim<3> >& derivs,
       const double multiplier,
       const double t,
       const double dt)
{
   // Get the mass density of the fluid.
   KeyType rhoKey(key.first, HydroFieldNames::massDensity);
   const Field<Dim<3>, Dim<3>::Scalar>& massDensity = 
      state.scalarField(rhoKey);

   // Get the electron temperature.  This assumes that Te = Ti = T.
   KeyType TeKey(key.first, HydroFieldNames::temperature);
   const Field<Dim<3>, Dim<3>::Scalar>& electronTemperature = 
      state.scalarField(TeKey);

   // Finally, get the resistivity.
   KeyType KKey(key.first, MHDFieldNames::resistivity);
   Field<Dim<3>, Dim<3>::Scalar>& resistivity = 
      state.scalarField(KKey);

   // Update the resistivity according to the Spitzer model.
   for (int i = 0; i < resistivity.nodeList().numInternalNodes(); ++i)
   {
      Dim<3>::Scalar Te = electronTemperature[i];
      Dim<3>::Scalar rho = massDensity[i];
      Dim<3>::Scalar R;
      if ((Te > 0.0) && (rho > mRhoMin))
      {
         R = std::max(mRmax, mC / std::pow(Te, 1.5));
      } // end if
      else
      {
         R = mRmax;
      } // end else

      resistivity[i] = R;
   } // end for

} // end update
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
bool
SpitzerResistivityUpdatePolicy::
operator==(const Spheral::UpdatePolicyBase<Dim<3>, Field<Dim<3>, Dim<3>::Scalar> >& rhs) const
{
   // We're only equal if the other guy is the same type.
   const SpitzerResistivityUpdatePolicy* rhsPtr = 
      dynamic_cast<const SpitzerResistivityUpdatePolicy*>(&rhs);
   if (rhsPtr == 0) {
      return false;
   } else {
      return true;
   }
}
//----------------------------------------------------------------------------

}
