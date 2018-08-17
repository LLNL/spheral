//---------------------------------Spheral++----------------------------------//
// StrengthModel -- The interface base class for strength models.
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrengthModel_hh__
#define __Spheral_StrengthModel_hh__

#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

// Forward declarations.
namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class StrengthModel {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;

  // Constructors, destructor.
  StrengthModel();
  virtual ~StrengthModel();

  //............................................................................
  // The generic interface we require all strength models to provide.
  virtual void shearModulus(FieldSpace::Field<Dimension, Scalar>& shearModulus,
                            const FieldSpace::Field<Dimension, Scalar>& density,
                            const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                            const FieldSpace::Field<Dimension, Scalar>& pressure) const = 0;

  virtual void yieldStrength(FieldSpace::Field<Dimension, Scalar>& yieldStrength,
                             const FieldSpace::Field<Dimension, Scalar>& density,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                             const FieldSpace::Field<Dimension, Scalar>& pressure,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrain,
                             const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate) const = 0;
  //............................................................................
  // Some strength models optionally provide the following methods.
  virtual bool providesSoundSpeed() const { return false; }
  virtual bool providesBulkModulus() const { return false; }
  virtual void soundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                          const FieldSpace::Field<Dimension, Scalar>& density,
                          const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                          const FieldSpace::Field<Dimension, Scalar>& pressure,
                          const FieldSpace::Field<Dimension, Scalar>& fluidSoundSpeed) const;

  virtual void bulkModulus(FieldSpace::Field<Dimension, Scalar>& bulkModulus,
                           const FieldSpace::Field<Dimension, Scalar>& massDensity,
                           const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void meltSpecificEnergy(FieldSpace::Field<Dimension, Scalar>& meltSpecificEnergy,
                                  const FieldSpace::Field<Dimension, Scalar>& density,
                                  const FieldSpace::Field<Dimension, Scalar>& specficThermalEnergy) const;

  virtual void coldSpecificEnergy(FieldSpace::Field<Dimension, Scalar>& coldSpecificEnergy,
                                  const FieldSpace::Field<Dimension, Scalar>& density,
                                  const FieldSpace::Field<Dimension, Scalar>& specficThermalEnergy) const;

protected:
  // The following individual methods are deprecated.
  virtual double shearModulus(const double density,
                              const double specificThermalEnergy,
                              const double pressure) const { VERIFY2(false, "Individual values for StrengthModel::shearModulus is deprecated."); }

  virtual double yieldStrength(const double density,
                               const double specificThermalEnergy,
                               const double pressure,
                               const double plasticStrain,
                               const double plasticStrainRate) const { VERIFY2(false, "Individual values for StrengthModel::yieldStrength is deprecated."); }

  virtual double soundSpeed(const double density,
                            const double specificThermalEnergy,
                            const double pressure,
                            const double fluidSoundSpeed) const { VERIFY2(false, "Individual values for StrengthModel::soundSpeed is deprecated."); }

private:
  //--------------------------- Private Interface ---------------------------//
  // No copying or assignment.
  StrengthModel(const StrengthModel&);
  StrengthModel& operator=(const StrengthModel&);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class StrengthModel;
  }
}

#endif

