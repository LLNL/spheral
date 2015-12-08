//---------------------------------Spheral++----------------------------------//
// IsothermalEquationOfState -- The isothermal equation of state:  P = cs^2 rho
//
// Created by JMO, Sat Nov  3 15:44:16 PDT 2007
//----------------------------------------------------------------------------//
#ifndef __Spheral_IsothermalEquationOfState_hh__
#define __Spheral_IsothermalEquationOfState_hh__

#include "EquationOfState.hh"

namespace Spheral {
namespace Material {

template<typename Dimension>
class IsothermalEquationOfState: public EquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  IsothermalEquationOfState(const double K,
                            const double mu,
                            const PhysicalConstants& constants,
                            const double minimumPressure,
                            const double maximumPressure,
                            const MaterialPressureMinType minPressureType);
  ~IsothermalEquationOfState();

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

  // Access the member data.
  double K() const;
  double molecularWeight() const;
  
  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void setExternalPressure(double P);

  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mK;
  double mCs;
  double mMolecularWeight;
  double mExternalPressure;

  // No default constructor, copying, or assignment.
  IsothermalEquationOfState();
  IsothermalEquationOfState(const IsothermalEquationOfState&);
  IsothermalEquationOfState& operator=(const IsothermalEquationOfState&);

  using EquationOfState<Dimension>::mConstants;
};

}
}

#ifndef __GCCXML__
#include "IsothermalEquationOfStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace Material {
    template<typename Dimension> class IsothermalEquationOfState;
  }
}

#endif
