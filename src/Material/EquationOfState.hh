//---------------------------------Spheral++----------------------------------//
// EquationOfState -- Abstract base class for the equation of state classes.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_EquationOfState_hh__
#define __Spheral_EquationOfState_hh__

#include <limits>
#include "PhysicalConstants.hh"
#include "Utilities/DBC.hh"

// Forward declarations.
namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace Material {

enum MaterialPressureMinType {
  PressureFloor = 0,
  ZeroPressure = 1,
};

template<typename Dimension>
class EquationOfState {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  EquationOfState(const PhysicalConstants& constants,
                  const double minimumPressure,
                  const double maximumPressure,
                  const MaterialPressureMinType minPressureType);

  virtual ~EquationOfState();

  // We require any equation of state to define the following methods for Fields.
  virtual void setPressure(FieldSpace::Field<Dimension, Scalar>& Pressure,
                           const FieldSpace::Field<Dimension, Scalar>& massDensity,
                           const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setTemperature(FieldSpace::Field<Dimension, Scalar>& temperature,
                              const FieldSpace::Field<Dimension, Scalar>& massDensity,
                              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setSpecificThermalEnergy(FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                                        const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                        const FieldSpace::Field<Dimension, Scalar>& temperature) const = 0;

  virtual void setSpecificHeat(FieldSpace::Field<Dimension, Scalar>& specificHeat,
                               const FieldSpace::Field<Dimension, Scalar>& massDensity,
                               const FieldSpace::Field<Dimension, Scalar>& temperature) const = 0;

  virtual void setSoundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                             const FieldSpace::Field<Dimension, Scalar>& massDensity,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setGammaField(FieldSpace::Field<Dimension, Scalar>& gamma,
			     const FieldSpace::Field<Dimension, Scalar>& massDensity,
			     const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setBulkModulus(FieldSpace::Field<Dimension, Scalar>& bulkModulus,
			     const FieldSpace::Field<Dimension, Scalar>& massDensity,
			     const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  // The set of constants defining our units.
  const PhysicalConstants& constants() const;

  // The min and max allowed pressures.
  double minimumPressure() const;
  double maximumPressure() const;
  void minimumPressure(const double x);
  void maximumPressure(const double x);
  
  // The algorithm for applying the minimum pressure.
  MaterialPressureMinType minimumPressureType() const;
  void minimumPressureType(const MaterialPressureMinType x);

  // Equations of state should have a valid test.
  virtual bool valid() const = 0;

  // Apply limits to a pressure value.
  Scalar applyPressureLimits(const Scalar P) const;

protected:
  // The following methods have been moved to the protected interface as we are deprecating them.
  virtual Scalar pressure(const Scalar massDensity,
                          const Scalar specificThermalEnergy) const { VERIFY2(false, "Individual values for EquationOfState::pressure is deprecated."); }

  virtual Scalar temperature(const Scalar massDensity,
                             const Scalar specificThermalEnergy) const { VERIFY2(false, "Individual values for EquationOfState::temperature is deprecated."); }

  virtual Scalar specificThermalEnergy(const Scalar massDensity,
                                       const Scalar temperature) const { VERIFY2(false, "Individual values for EquationOfState::specificThermalEnergy is deprecated."); }

  virtual Scalar specificHeat(const Scalar massDensity,
                              const Scalar temperature) const { VERIFY2(false, "Individual values for EquationOfState::specificHeat is deprecated."); }

  virtual Scalar soundSpeed(const Scalar massDensity,
                            const Scalar specificThermalEnergy) const { VERIFY2(false, "Individual values for EquationOfState::soundSpeed is deprecated."); }

  virtual Scalar gamma(const Scalar massDensity,
		       const Scalar specificThermalEnergy) const { VERIFY2(false, "Individual values for EquationOfState::gamma is deprecated."); }

  virtual Scalar bulkModulus(const Scalar massDensity,
                             const Scalar specificThermalEnergy) const { VERIFY2(false, "Individual values for EquationOfState::bulkModulus is deprecated."); }

  PhysicalConstants mConstants;

private:
  double mMinimumPressure, mMaximumPressure;
  MaterialPressureMinType mMinPressureType;

  // No default constructor.
  EquationOfState();
};

}
}

#ifndef __GCCXML__
#include "EquationOfStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace Material {
    template<typename Dimension> class EquationOfState;
  }
}

#endif
