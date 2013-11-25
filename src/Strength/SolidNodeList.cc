//---------------------------------Spheral++----------------------------------//
// SolidNodeList -- A form of the SPH NodeList appropriate for use with 
// solid materials.
//
// Created by JMO, Tue Sep 7 22:44:37 2004
//----------------------------------------------------------------------------//
#include "SolidNodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Material/EquationOfState.hh"
#include "Kernel/TableKernel.hh"
#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "SolidMaterial/StrengthModel.hh"

#include "Utilities/DBC.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;

using Material::EquationOfState;
using FieldSpace::Field;
using FieldSpace::FieldList;
using FileIOSpace::FileIO;
using NodeSpace::FluidNodeList;

//------------------------------------------------------------------------------
// Construct with the given EOS object, along with optional numInternal nodes,
// numGhost nodes, and name.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidNodeList<Dimension>::
SolidNodeList(string name,
              EquationOfState<Dimension>& eos,
              SolidMaterial::StrengthModel<Dimension>& strength,
              const int numInternal,
              const int numGhost,
              const Scalar hmin,
              const Scalar hmax,
              const Scalar hminratio,
              const Scalar nPerh,
              const int maxNumNeighbors,
              const Scalar rhoMin,
              const Scalar rhoMax):
  FluidNodeList<Dimension>(name, 
                           eos,
                           numInternal, 
                           numGhost,
                           hmin,
                           hmax,
                           hminratio,
                           nPerh,
                           maxNumNeighbors,
                           rhoMin,
                           rhoMax),
  mDeviatoricStress(SolidFieldNames::deviatoricStress, *this),
  mPlasticStrain(SolidFieldNames::plasticStrain, *this),
  mPlasticStrainRate(SolidFieldNames::plasticStrainRate, *this),
  mDamage(SolidFieldNames::tensorDamage, *this),
  mEffectiveDamage(SolidFieldNames::effectiveTensorDamage, *this),
  mDamageGradient(SolidFieldNames::damageGradient, *this),
  mFragmentIDs(SolidFieldNames::fragmentIDs, *this),
  mStrength(strength) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidNodeList<Dimension>::
~SolidNodeList() {
}

//------------------------------------------------------------------------------
// Calculate and return the sound speed field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidNodeList<Dimension>::
soundSpeed(Field<Dimension, typename Dimension::Scalar>& field) const {

  // Get the straight EOS look up from the base class.
  FluidNodeList<Dimension>::soundSpeed(field);

  // Augment the sound speed with the strength model.
  const Field<Dimension, Scalar>& rho = this->massDensity();
  const Field<Dimension, Scalar>& u = this->specificThermalEnergy();
  Field<Dimension, Scalar> P(HydroFieldNames::pressure, *this);
  this->pressure(P);
  for (int i = 0; i != this->numInternalNodes(); ++i) {
    field(i) = mStrength.soundSpeed(rho(i), u(i), P(i), field(i));
  }
}

//------------------------------------------------------------------------------
// Calculate and return the bulk modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidNodeList<Dimension>::
bulkModulus(Field<Dimension, typename Dimension::Scalar>& field) const {
  const Field<Dimension, Scalar>& rho = this->massDensity();
  const Field<Dimension, Scalar>& u = this->specificThermalEnergy();
  const EquationOfState<Dimension>& eos = this->equationOfState();
  this->equationOfState().setBulkModulus(field, rho, u);
}

//------------------------------------------------------------------------------
// Calculate and return the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidNodeList<Dimension>::
shearModulus(Field<Dimension, typename Dimension::Scalar>& field) const {
  const Field<Dimension, Scalar>& rho = this->massDensity();
  const Field<Dimension, Scalar>& u = this->specificThermalEnergy();
  Field<Dimension, Scalar> P(HydroFieldNames::pressure, *this);
  this->pressure(P);
  mStrength.shearModulus(field, rho, u, P);
}

//------------------------------------------------------------------------------
// Calculate and return the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidNodeList<Dimension>::
yieldStrength(Field<Dimension, typename Dimension::Scalar>& field) const {
  const Field<Dimension, Scalar>& rho = this->massDensity();
  const Field<Dimension, Scalar>& u = this->specificThermalEnergy();
  Field<Dimension, Scalar> P(HydroFieldNames::pressure, *this);
  this->pressure(P);
  const Field<Dimension, SymTensor>& D = this->effectiveDamage();
  mStrength.yieldStrength(field, rho, u, P, mPlasticStrain, mPlasticStrainRate);
}

//------------------------------------------------------------------------------
// Dump the current state of the NodeList to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidNodeList<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Dump the ancestor class.
  FluidNodeList<Dimension>::dumpState(file, pathName);

  file.write(mDeviatoricStress, pathName + "/" + mDeviatoricStress.name());
  file.write(mPlasticStrain, pathName + "/" + mPlasticStrain.name());
  file.write(mPlasticStrainRate, pathName + "/" + mPlasticStrainRate.name());
  file.write(mDamage, pathName + "/" + mDamage.name());
  file.write(mEffectiveDamage, pathName + "/" + mEffectiveDamage.name());
  file.write(mDamageGradient, pathName + "/" + mDamageGradient.name());
  file.write(mFragmentIDs, pathName + "/" + mFragmentIDs.name());
}

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidNodeList<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // Restore the ancestor class.
  FluidNodeList<Dimension>::restoreState(file, pathName);

  file.read(mDeviatoricStress, pathName + "/" + mDeviatoricStress.name());
  file.read(mPlasticStrain, pathName + "/" + mPlasticStrain.name());
  file.read(mPlasticStrainRate, pathName + "/" + mPlasticStrainRate.name());
  file.read(mDamage, pathName + "/" + mDamage.name());
  file.read(mEffectiveDamage, pathName + "/" + mEffectiveDamage.name());
  file.read(mDamageGradient, pathName + "/" + mDamageGradient.name());
  file.read(mFragmentIDs, pathName + "/" + mFragmentIDs.name());
}

}
}

