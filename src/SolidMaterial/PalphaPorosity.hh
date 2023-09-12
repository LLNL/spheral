//---------------------------------Spheral++----------------------------------//
// PalphaPorosity.hh
//
// An implementation of the P-alpha porosity model described in
//
// Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
// involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
// hydrocode. Icarus, 198(1), 242–255.
//
// This model assumes you will provide a solid EOS which will be modified.
// The underlying actualy solid EOS should provide the reference density, which
// will be treated here as the compacted true solid reference density.
//
// Note this model introduces a new state variable, the distention (alpha), which
// the pressure now depends on.  This implies our usual definition of P(rho, eps)
// now becomes P(rho, eps, alpha).  Our EOS interface does not recognize this
// this parameter, so we store alpha locally and only allow Field updates of the
// pressure (forbidding the single value P lookup the EOS usually allows).
//
// PalphaPorosity is the physics module which time evolves the distention 
// parameter (alpha) and gives it to the PorousEquationOfState.
//
// Created by JMO, Wed Sep  6 10:08:15 PDT 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PalphaPorosity__
#define __Spheral_PalphaPorosity__

#include "PorousEquationOfState.hh"
#include "PorousStrengthModel.hh"
#include "Physics/Physics.hh"
#include "DataOutput/registerWithRestart.hh"

#include <limits>

namespace Spheral {

template<typename Dimension>
class PalphaPorosity: 
    public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors, destructors.
  PalphaPorosity(PorousEquationOfState<Dimension>& porousEOS,     // Porous EOS we're going to modify
                 PorousStrengthModel<Dimension>& porousStrength,  // Porous strength model we're going to modify
                 const NodeList<Dimension>& nodeList,             // The NodeList we're going apply to
                 const double phi0,                               // Initial porosity
                 const double Pe,                                 // Elastic pressure threshold
                 const double Pt,                                 // Transition pressure (Pe <= Pt)
                 const double alphae,                             // Elastic distension threshold
                 const double alphat,                             // Transition distension
                 const double n1,                                 // Fitted exponent for plastic distention evolution
                 const double n2,                                 // Fitted exponent for plastic distention evolution
                 const double cS0,                                // Reference sound speed at full density
                 const double c0);                                // Reference sound speed at initial porosity

  PalphaPorosity(PorousEquationOfState<Dimension>& porousEOS,     // Porous EOS we're going to modify
                 PorousStrengthModel<Dimension>& porousStrength,  // Porous strength model we're going to modify
                 const NodeList<Dimension>& nodeList,             // The NodeList we're going apply to
                 const Field<Dimension, Scalar>& phi0,            // Initial porosity
                 const double Pe,                                 // Elastic pressure threshold
                 const double Pt,                                 // Transition pressure (Pe <= Pt)
                 const double alphae,                             // Elastic distension threshold
                 const double alphat,                             // Transition distension
                 const double n1,                                 // Fitted exponent for plastic distention evolution
                 const double n2,                                 // Fitted exponent for plastic distention evolution
                 const double cS0,                                // Reference sound speed at full density
                 const Field<Dimension, Scalar>& c0);             // Reference sound speed at initial porosity

  virtual ~PalphaPorosity();

  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);

  // Do any required one-time initializations on problem start up.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase);
  //............................................................................

  //............................................................................
  // Methods required for restarting.
  virtual std::string label() const { return "PalphaPorosity " + mNodeList.name(); }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //............................................................................

  // Access the material parameters.
  double Pe() const;
  double Pt() const;
  double alphae() const;
  double alphat() const;
  double n1() const;
  double n2() const;
  double cS0() const;
  const PorousEquationOfState<Dimension>& porousEOS() const;
  const PorousStrengthModel<Dimension>& porousStrength() const;
  const NodeList<Dimension>& nodeList() const;
  const Field<Dimension, Scalar>& c0() const;
  const Field<Dimension, Scalar>& alpha0() const;
  const Field<Dimension, Scalar>& alpha() const;
  const Field<Dimension, Scalar>& DalphaDt() const;
  const Field<Dimension, Scalar>& partialPpartialEps() const;
  const Field<Dimension, Scalar>& partialPpartialRho() const;

  // Provide the porosity (phi) computed from the internally stored distention alpha
  Field<Dimension, Scalar> phi() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mPe, mPt, mAlphae, mAlphat, mn1, mn2, mcS0;
  PorousEquationOfState<Dimension>& mPorousEOS;
  PorousStrengthModel<Dimension>& mPorousStrength;
  const NodeList<Dimension>& mNodeList;
  Field<Dimension, Scalar> mc0, mAlpha0, mAlpha, mDalphaDt, mdPdU, mdPdR;

  // The restart registration.
  RestartRegistrationType mRestart;

  // Disallow default constructor
  PalphaPorosity();
};

}

#include "PalphaPorosityInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PalphaPorosity;
}

#endif
