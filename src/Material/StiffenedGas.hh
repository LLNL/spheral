//---------------------------------Spheral++----------------------------------//
// StiffenedGas -- stiffened gamma gas law
//
//----------------------------------------------------------------------------//

#ifndef StiffenedGas_HH
#define StiffenedGas_HH

#include "EquationOfState.hh"

namespace Spheral {

template<typename Dimension>
class StiffenedGas: public EquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  StiffenedGas(const double gamma,
               const double P0, 
               const double Cv,
               const PhysicalConstants& constants,
               const double minimumPressure,
               const double maximumPressure,
               const MaterialPressureMinType minPressureType,
               const double externalPressure);
  ~StiffenedGas();

  // We require any equation of state to define the following properties.
  virtual void setPressure(Field<Dimension, Scalar>& Pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const override;

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const override;

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setGammaField(Field<Dimension, Scalar>& gamma,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override;

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const override;

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
  Scalar gamma() const;
  void gamma(Scalar gamma);

  Scalar specificHeat() const;
  void specificHeat(Scalar Cv);

  Scalar referencePressure() const;
  void referencePressure(Scalar P0);
  
  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
  double mGamma;
  double mGamma1;
  double mP0;
  double mCv;

  // No default constructor, copying, or assignment.
  StiffenedGas();
  StiffenedGas(const StiffenedGas&);
  StiffenedGas& operator=(const StiffenedGas&);
  
  using EquationOfState<Dimension>::mConstants;
};

}

#include "StiffenedGasInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class StiffenedGas;
}

#endif
