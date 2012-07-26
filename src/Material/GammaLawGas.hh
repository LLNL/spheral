//---------------------------------Spheral++----------------------------------//
// GammaLawGas -- The gamma law gas equation of state.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//

#ifndef GammaLawGas_HH
#define GammaLawGas_HH

#include "EquationOfState.hh"

namespace Spheral {
namespace Material {

template<typename Dimension, typename Constants>
class GammaLawGas: public EquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  GammaLawGas(const double gamma,
              const double mu,
              const double minimumPressure = -std::numeric_limits<double>::max(),
              const double maximumPressure = std::numeric_limits<double>::max());
  ~GammaLawGas();

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
  double getGamma() const;
  void setGamma(double gamma);

  double getMolecularWeight() const;
  void setMolecularWeight(double molecularWeight);
  
  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mGamma;
  double mGamma1;
  double mMolecularWeight;

  // No default constructor, copying, or assignment.
  GammaLawGas();
  GammaLawGas(const GammaLawGas&);
  GammaLawGas& operator=(const GammaLawGas&);
  
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace Material {
    template<typename Dimension, typename Constants> class GammaLawGas;
  }
}

#endif
