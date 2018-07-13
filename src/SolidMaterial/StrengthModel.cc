//---------------------------------Spheral++----------------------------------//
// StrengthModel -- The interface base class for strength models.
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//

#include "StrengthModel.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace SolidMaterial {

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
StrengthModel<Dimension>::
StrengthModel() {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
StrengthModel<Dimension>::
~StrengthModel() {
}

//------------------------------------------------------------------------------
// Melt energy
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrengthModel<Dimension>::
meltSpecificThermalEnergy(FieldSpace::Field<Dimension, Scalar>& meltSpecificEnergy,
                          const FieldSpace::Field<Dimension, Scalar>& density,
                          const FieldSpace::Field<Dimension, Scalar>& specficThermalEnergy) const {
  meltSpecificEnergy = 0.0;
}

//------------------------------------------------------------------------------
// Cold energy
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrengthModel<Dimension>::
coldSpecificThermalEnergy(FieldSpace::Field<Dimension, Scalar>& coldSpecificEnergy,
                          const FieldSpace::Field<Dimension, Scalar>& density,
                          const FieldSpace::Field<Dimension, Scalar>& specficThermalEnergy) const {
  coldSpecificEnergy = 0.0;
}

}
}
