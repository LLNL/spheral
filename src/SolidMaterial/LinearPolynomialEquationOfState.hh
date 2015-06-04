//---------------------------------Spheral++----------------------------------//
// LinearPolynomialEquationOfState -- An equation of state approximated by a
// linear polynomial, i.e.:
//
//   P(rho, e) = A0 + A1*mu + a2*mu^2 + a3*mu^3 + (B0 + B1*mu + B2*mu^2)*e
//   mu = rho/rho0 - 1.0
//
// Created by JMO, Thu May  5 16:07:36 PDT 2005
//----------------------------------------------------------------------------//
#ifndef LinearPolynomialEquationOfState_HH
#define LinearPolynomialEquationOfState_HH

#include <float.h>

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
class LinearPolynomialEquationOfState: public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  LinearPolynomialEquationOfState(const double referenceDensity,
                                  const double etamin,
                                  const double etamax,
                                  const double a0,
                                  const double a1,
                                  const double a2,
                                  const double a3,
                                  const double b0,
                                  const double b1,
                                  const double b2,
                                  const double atomicWeight,
                                  const Material::PhysicalConstants& constants,
                                  const double externalPressure,
                                  const double minimumPressure,
                                  const double maximumPressure,
                                  const Material::MaterialPressureMinType minPressureType);
  ~LinearPolynomialEquationOfState();

  // We require any equation of state to define the following properties.
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
  double a0() const;
  double a1() const;
  double a2() const;
  double a3() const;
  double b0() const;
  double b1() const;
  double b2() const;
  double atomicWeight() const;
  
  void a0(const double x);
  void a1(const double x);
  void a2(const double x);
  void a3(const double x);
  void b0(const double x);
  void b1(const double x);
  void b2(const double x);
  void atomicWeight(const double x);
  
  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void externalPressure(const double x);

  // Compute the derivative of the pressure with respect to the density.
  double computeDPDrho(const Scalar massDensity,
                       const Scalar specificThermalEnergy) const;

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

  virtual Scalar gamma(const Scalar massDensity,
		       const Scalar specificThermalEnergy) const;

  virtual Scalar bulkModulus(const Scalar massDensity,
                             const Scalar specificThermalEnergy) const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mA0, mA1, mA2, mA3;
  double mB0, mB1, mB2;
  double mAtomicWeight;
  double mCv;
  double mGamma;
  double mExternalPressure;

  // Disallow default constructor
  LinearPolynomialEquationOfState();

};

}
}

#ifndef __GCCXML__
#include "LinearPolynomialEquationOfStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class LinearPolynomialEquationOfState;
  }
}

#endif
