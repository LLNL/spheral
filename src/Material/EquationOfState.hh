//---------------------------------Spheral++----------------------------------//
// EquationOfState -- Abstract base class for the equation of state classes.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_EquationOfState_hh__
#define __Spheral_EquationOfState_hh__

#include "PhysicalConstants.hh"
#include "Utilities/DBC.hh"

#include <limits>

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class Field;

enum class MaterialPressureMinType {
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
                  const MaterialPressureMinType minPressureType,
                  const double externalPressure);

  virtual ~EquationOfState();

  // We require any equation of state to define the following methods for Fields.
  virtual void setPressure(FieldView<Dimension, Scalar>& Pressure,
                           const FieldView<Dimension, Scalar>& massDensity,
                           const FieldView<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setPressureAndDerivs(FieldView<Dimension, Scalar>& Pressure,           // set pressure
                                    FieldView<Dimension, Scalar>& dPdu,               // set (\partial P)/(\partial u) (specific thermal energy)
                                    FieldView<Dimension, Scalar>& dPdrho,             // set (\partial P)/(\partial rho) (density)
                                    const FieldView<Dimension, Scalar>& massDensity,
                                    const FieldView<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const = 0;

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const = 0;

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setGammaField(Field<Dimension, Scalar>& gamma,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const = 0;

  // Look up an energy that gives the requested pressure at the specified density.
  virtual Scalar specificThermalEnergyForPressure(const Scalar Ptarget,
                                                  const Scalar rho,
                                                  const Scalar epsMin,
                                                  const Scalar epsMax,
                                                  const Scalar epsTol,
                                                  const Scalar Ptol,
                                                  const unsigned maxIterations,
                                                  const bool verbose = false) const;

  // Optionally provide a molecular weight.
  virtual Scalar molecularWeight() const;

  // The set of constants defining our units.
  const PhysicalConstants& constants() const;

  // The min and max allowed pressures.
  double minimumPressure() const;
  double maximumPressure() const;
  void minimumPressure(double x);
  void maximumPressure(double x);
  
  // The algorithm for applying the minimum pressure.
  MaterialPressureMinType minimumPressureType() const;
  void minimumPressureType(MaterialPressureMinType x);

  // Optionally we can provide an "external" pressure, which is just a pressure offset
  double externalPressure() const;
  void externalPressure(const double x);

  // Equations of state should have a valid test.
  virtual bool valid() const = 0;

  // Apply limits to a pressure value.
  Scalar applyPressureLimits(Scalar P) const;

protected:
  PhysicalConstants mConstants;

private:
  double mMinimumPressure, mMaximumPressure, mExternalPressure;
  MaterialPressureMinType mMinPressureType;

  // No default constructor.
  EquationOfState();
};

}

#include "EquationOfStateInline.hh"

#endif
