//---------------------------------Spheral++----------------------------------//
// GruneisenEquationOfState -- Gruneisen  equation of state.
//
// Created by MLF & JMO Wed Sep 17 13:32:54 EDT 2003 
// Reference: Equation of State and Strength of Properties of Selected Materials
//            Daniel J. Steinberg, UCRL-MA-106439, February 13, 1991
//----------------------------------------------------------------------------//
#ifndef GruneisenEquationOfState_HH
#define GruneisenEquationOfState_HH

#include <limits>
#include "SolidEquationOfState.hh"

// Forward declarations.
namespace Spheral {
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class GruneisenEquationOfState: public SolidEquationOfState<Dimension> { 

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  GruneisenEquationOfState(const double referenceDensity,
                           const double etamin,
                           const double etamax,
                           const double C0, 
                           const double S1,
                           const double S2,
                           const double S3,
                           const double gamma0,
                           const double b,
                           const double atomicWeight,
                           const Material::PhysicalConstants& constants,
                           const double externalPressure,
                           const double minimumPressure,
                           const double maximumPressure,
                           const Material::MaterialPressureMinType minPressureType);
  virtual ~GruneisenEquationOfState();

  // We require any equation of state to define the following methods for Fields.
  virtual void setPressure(FieldSpace::Field<Dimension, Scalar>& Pressure,
                           const FieldSpace::Field<Dimension, Scalar>& massDensity,
                           const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setTemperature(FieldSpace::Field<Dimension, Scalar>& temperature,
                              const FieldSpace::Field<Dimension, Scalar>& massDensity,
                              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setSpecificThermalEnergy(FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
                                        const FieldSpace::Field<Dimension, Scalar>& massDensity,
                                        const FieldSpace::Field<Dimension, Scalar>& temperature) const;

  virtual void setSpecificHeat(FieldSpace::Field<Dimension, Scalar>& specificHeat,
                               const FieldSpace::Field<Dimension, Scalar>& massDensity,
                               const FieldSpace::Field<Dimension, Scalar>& temperature) const;

  virtual void setSoundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
                             const FieldSpace::Field<Dimension, Scalar>& massDensity,
                             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setGammaField(FieldSpace::Field<Dimension, Scalar>& gamma,
			     const FieldSpace::Field<Dimension, Scalar>& massDensity,
			     const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  virtual void setBulkModulus(FieldSpace::Field<Dimension, Scalar>& bulkModulus,
                              const FieldSpace::Field<Dimension, Scalar>& massDensity,
                              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  // Access the member data.
  double C0() const;
  double S1() const;
  double S2() const;
  double S3() const;
  double gamma0() const;
  double b() const;
  double atomicWeight() const;
  double Cv() const;

  void C0(const double val);
  void S1(const double val);
  void S2(const double val);
  void S3(const double val);
  void gamma0(const double val);
  void b(const double val);
  void atomicWeight(const double val);

  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void externalPressure(const double P);

  // Compute (\partial P)/(\partial rho) for use in sound speed and bulk modulus.
  Scalar computeDPDrho(const Scalar massDensity,
                       const Scalar specificThermalEnergy) const;

  // Equations of state should have a valid test.
  virtual bool valid() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // We also want the equivalent functions for individual calculations.
  virtual Scalar pressure(const Scalar massDensity,
                          const Scalar specificThermalEnergy) const;

  virtual Scalar temperature(const Scalar massDensity,
                             const Scalar specificThermalEnergy) const;

  virtual Scalar specificThermalEnergy(const Scalar massDensity,
                                       const Scalar temperature) const;

  virtual Scalar specificHeat(const Scalar massDensity,
                              const Scalar temperature) const;

  virtual Scalar soundSpeed(const Scalar massDensity,
                            const Scalar specificThermalEnergy) const;

  virtual Scalar gamma(const Scalar massdensity,
                       const Scalar specificThermalEnergy) const;

  virtual Scalar bulkModulus(const Scalar massDensity,
                             const Scalar specificThermalEnergy) const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mC0;
  double mS1;
  double mS2;
  double mS3;
  double mgamma0;
  double mb;
  double mAtomicWeight;
  double mCv;
  double mExternalPressure;

  // No default constructor, copying, or assignment.
  GruneisenEquationOfState();
  GruneisenEquationOfState(const GruneisenEquationOfState&);
  GruneisenEquationOfState& operator=(const GruneisenEquationOfState&);

  using Material::EquationOfState<Dimension>::mConstants;
};

}
}

#ifndef __GCCXML__
#include "GruneisenEquationOfStateInline.hh"
#endif

#else
// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class GruneisenEquationOfState;
  }
}

#endif
