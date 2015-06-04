//---------------------------------Spheral++----------------------------------//
// TillotsonEquationOfState -- This equation of state is designed to represent
// metallic materials over a range of pressure and density -- spanning solid, 
// liquid, and vapor states.
//
// Tillotson 1962
//
// Created by JMO, Wed Mar 16 23:31:17 PDT 2011
//----------------------------------------------------------------------------//
#ifndef TillotsonEquationOfState_HH
#define TillotsonEquationOfState_HH

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
class TillotsonEquationOfState: public SolidEquationOfState<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructors.
  TillotsonEquationOfState(const double referenceDensity,
                           const double etamin,
                           const double etamax,
                           const double etamin_solid,
                           const double etamax_solid,
                           const double a,
                           const double b,
                           const double A,
                           const double B,
                           const double alpha,
                           const double beta,
                           const double eps0,
                           const double epsLiquid,
                           const double epsVapor,
                           const double atomicWeight,
                           const Material::PhysicalConstants& constants,
                           const double externalPressure,
                           const double minimumPressure,
                           const double maximumPressure,
                           const Material::MaterialPressureMinType minPressureType);
  ~TillotsonEquationOfState();

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
  double etamin_solid() const;
  double etamax_solid() const;
  double a() const;
  double b() const;
  double A() const;
  double B() const;
  double alpha() const;
  double beta() const;
  double eps0() const;
  double epsLiquid() const;
  double epsVapor() const;
  double atomicWeight() const;
  
  void etamin_solid(const double x);
  void etamax_solid(const double x);
  void a(const double x);
  void b(const double x);
  void A(const double x);
  void B(const double x);
  void alpha(const double x);
  void beta(const double x);
  void eps0(const double x);
  void epsLiquid(const double x);
  void epsVapor(const double x);
  void atomicWeight(const double x);
  
  // If requested, the user can specify an external pressure to be applied
  // in the pressure calculation.
  double externalPressure() const;
  void externalPressure(const double x);

  // Compute the derivative of the pressure with respect to the density.
  double computeDPDrho(const Scalar massDensity,
                       const Scalar specificThermalEnergy) const;

  // Helper methods for computing the components of the pressure and it's derivatives.
  double computePhi(const double& eta, const double& eps) const;
  double computeP1(const double& mu, const double& P2) const;
  double computeP2(const double& phi, const double& mu, const double& rho, const double& eps) const;
  double computeP4(const double& phi, const double& mu, const double& eta, const double& rho, const double& eps) const;
  double compute_dphidrho_eps(const double& rho0, const double& eta, const double& eps) const;
  double compute_dP1drho_eps(const double& rho0, const double& mu, const double& dP2drho_eps) const;
  double compute_dP2drho_eps(const double& phi, const double& dphidrho_eps, const double& rho0, const double& rho, const double& eps) const;
  double compute_dP4drho_eps(const double& phi, const double& dphidrho_eps, const double& rho0, const double& eta, const double& mu, const double& rho, const double& eps) const;
  double compute_dphideps_rho(const double& eta, const double& eps) const;
  double compute_dP2deps_rho(const double& phi, const double& dphideps_rho, const double& rho, const double& eps) const;
  double compute_dP4deps_rho(const double& phi, const double& dphideps_rho, const double& eta, const double& rho, const double& eps) const;

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
  double mEtaMinSolid, mEtaMaxSolid,
         ma, mb, mA, mB, malpha, mbeta, meps0, mepsLiquid, mepsVapor,
         mAtomicWeight, mCv, mGamma, mExternalPressure;

  // Disallow default constructor
  TillotsonEquationOfState();

  using Material::EquationOfState<Dimension>::mConstants;
};

}
}

#ifndef __GCCXML__
#include "TillotsonEquationOfStateInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class TillotsonEquationOfState;
  }
}

#endif
