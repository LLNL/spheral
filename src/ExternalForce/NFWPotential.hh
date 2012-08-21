//---------------------------------Spheral++----------------------------------//
// NFWPotential -- Impose a potential due to a Navarro Frenk White mass 
// density profile, which is supposed to be a "universal" dark matter halo
// profile seen in cosmological N body simulations.
// Navarro, Frenk, & White 1997, ApJ, 490, 493-508.
//
// Note that this object is templated on your unit choice, since we depend on
// G and time units (via how the Hubble constant is defined.)
//
// NFW density profile:
// rho/rho_0(r) = delta_c/{ (r/r_s)*(1 + r/r_s)^2 }
//
// Corresponding mass profile:
// M(r) = 4/3 pi delta_c rho_0 r_s^4 { 1 + r/r_s - 2 ln(1 + r/r_s) - 1/(1 + r/r_s) }
//
// Balancing velocity profile:
// v^2(r) = 4/3 pi delta_c rho_0 r_s^4 G {( 1/(1 + r/r_s) + 2 ln(1 + r/r_s) - 1 )/r +
//                                        1/(1 + r/r_s)^2 - 2/(1 + r/r_s) }
//
// Created by JMO, Thu May  8 18:02:30 PDT 2003
//----------------------------------------------------------------------------//
#ifndef NFWPotential_HH
#define NFWPotential_HH

#include "Physics/GenericBodyForce.hh"
#include "Material/PhysicalConstants.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace Material {
    class PhysicalConstants;
  }
}

namespace Spheral {
namespace PhysicsSpace {

template<typename Dimension>
class NFWPotential: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  // deltac = characteristic dimensionless density
  // rs = scale radius
  // h0 = Hubble constant (units of 100 km/sec/Mpc)
  NFWPotential(double deltac,
               double rs,
               double h0,
               const Vector& origin,
               const Material::PhysicalConstants& constants);

  // Destructor.
  virtual ~NFWPotential();

  // This is the derivative method that all BodyPotential classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Provide the timestep appropriate for this package.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Get the cumulative potential energy calculated in the last 
  // evaluateDerivatives.
  virtual Scalar extraEnergy() const;

  //! Required label for Physics interface.
  virtual std::string label() const { return "NFWPotential"; }

  // Calculate the mass density as a function of radius.
  Scalar massDensity(Scalar r) const;

  // Calculate the enclosed mass as a function of radius.
  Scalar enclosedMass(Scalar r) const;

  // Calculate the supporting orbital velocity as a function of radius.
  Scalar orbitalVelocity(Scalar r) const;

  // The characteristic density.
  Scalar characteristicDensity() const;
  void setCharacteristicDensity(const Scalar x);

  // The scale Radius.
  Scalar scaleRadius() const;
  void setScaleRadius(const Scalar x);

  // The hubble constant (in units of 100 km/sec/Mpc).
  Scalar h0() const;
  void seth0(const Scalar rc);

  // The critical (closure) density for the universe.
  Scalar criticalDensity() const;

  // The origin (center of the profile).
  const Vector& origin() const;
  void setOrigin(const Vector& origin);

  // The maximum allowed fractional change in a particles potential (for 
  // setting the timestep).
  Scalar deltaPotentialFraction() const;
  void setDeltaPotentialFraction(const Scalar deltaPhi);

private:
  //--------------------------- Public Interface ---------------------------//
  Scalar mDeltac;
  Scalar mRs;
  Scalar mh0;
  Vector mOrigin;
  Material::PhysicalConstants mConstants;
  Scalar mDeltaPhiFraction;
  Scalar mCriticalDensity;
  mutable Scalar mPotentialEnergy;


  // No default constructor, copying, or assignment.
  NFWPotential();
  NFWPotential(const NFWPotential& rhs);
  NFWPotential& operator=(const NFWPotential& rhs);
};

}
}

#ifndef __GCCXML__
#include "NFWPotentialInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class NFWPotential;
  }
}

#endif
