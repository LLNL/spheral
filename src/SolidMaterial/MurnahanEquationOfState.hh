//---------------------------------Spheral++----------------------------------//
// MurnahanEquationOfState
//
//   P(rho) = 1/(nK) * (eta^n - 1)
//   eta = rho/rho0
//
// Created by JMO, Mon Jun  6 13:53:50 PDT 2005
//----------------------------------------------------------------------------//
#ifndef MurnahanEquationOfState_HH
#define MurnahanEquationOfState_HH

#include <float.h>
#include "SolidEquationOfState.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class MurnahanEquationOfState: public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  MurnahanEquationOfState(const double referenceDensity,
                          const double etamin,
                          const double etamax,
                          const double n,
                          const double K,
                          const double atomicWeight,
                          const Material::PhysicalConstants& constants,
                          const double externalPressure,
                          const double minimumPressure,
                          const double maximumPressure,
                          const Material::MaterialPressureMinType minPressureType);
  ~MurnahanEquationOfState();

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

  virtual void setEntropy(FieldSpace::Field<Dimension, Scalar>& entropy,
                          const FieldSpace::Field<Dimension, Scalar>& massDensity,
                          const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const;

  // We also want the equivalent functions for individual calculations.
  Scalar pressure(const Scalar massDensity,
                  const Scalar specificThermalEnergy) const;

  Scalar temperature(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar specificThermalEnergy(const Scalar massDensity,
                               const Scalar temperature) const;

  Scalar specificHeat(const Scalar massDensity,
                      const Scalar temperature) const;

  Scalar soundSpeed(const Scalar massDensity,
                    const Scalar specificThermalEnergy) const;

  Scalar gamma(const Scalar massDensity,
               const Scalar specificThermalEnergy) const;

  Scalar bulkModulus(const Scalar massDensity,
                     const Scalar specificThermalEnergy) const;

  Scalar entropy(const Scalar massDensity,
                 const Scalar specificThermalEnergy) const;

  // Access the member data.
  double n() const;
  double K() const;
  double atomicWeight() const;
  
  void n(const double x);
  void K(const double x);
  void atomicWeight(const double x);
  
  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void externalPressure(const double x);

  // Compute the derivative of the pressure with respect to the density.
  double computeDPDrho(const Scalar massDensity,
                       const Scalar specificThermalEnergy) const;

  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mn, mK;
  double mAtomicWeight;
  double mExternalPressure;
  double mCv;
  double mnKi;

  // Disallow default constructor
  MurnahanEquationOfState();

  using Material::EquationOfState<Dimension>::mConstants;
};

}
}

#ifndef __GCCXML__
#include "MurnahanEquationOfStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class MurnahanEquationOfState;
  }
}

#endif
