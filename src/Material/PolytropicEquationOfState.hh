//---------------------------------Spheral++----------------------------------//
// PolytropicEquationOfState -- A polytrope equation of state.
//
// Created by JMO, Wed Apr  2 14:48:12 PST 2003
//----------------------------------------------------------------------------//
#ifndef __Spheral_PolytropicEquationOfState_hh__
#define __Spheral_PolytropicEquationOfState_hh__

#include "EquationOfState.hh"

namespace Spheral {
namespace Material {

template<typename Dimension>
class PolytropicEquationOfState: public EquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  PolytropicEquationOfState(const double K,
                            const double index,
                            const double mu,
                            const PhysicalConstants& constants,
                            const double minimumPressure,
                            const double maximumPressure,
                            const MaterialPressureMinType minPressureType);
  ~PolytropicEquationOfState();

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

  // Access the member data.
  double polytropicConstant() const;
  double polytropicIndex() const;
  double gamma() const;
  double molecularWeight() const;
  
  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void setExternalPressure(double P);

  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mPolytropicConstant;
  double mPolytropicIndex;
  double mGamma;
  double mGamma1;
  double mMolecularWeight;

  double mExternalPressure;

  // No default constructor, copying, or assignment.
  PolytropicEquationOfState();
  PolytropicEquationOfState(const PolytropicEquationOfState&);
  PolytropicEquationOfState& operator=(const PolytropicEquationOfState&);

  using EquationOfState<Dimension>::mConstants;
};

}
}

#ifndef __GCCXML__
#include "PolytropicEquationOfStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace Material {
    template<typename Dimension> class PolytropicEquationOfState;
  }
}

#endif
