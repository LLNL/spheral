//---------------------------------Spheral++----------------------------------//
// PalphaPorosity.hh
//
// An implementation of the P-alpha porosity model described in
//
// Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
// involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
// hydrocode. Icarus, 198(1), 242–255.
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

#include "Porosity/PorosityModel.hh"
#include "DataOutput/registerWithRestart.hh"

namespace Spheral {

template<typename Dimension>
class PalphaPorosity: 
    public PorosityModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors, destructors.
  PalphaPorosity(const SolidNodeList<Dimension>& nodeList,        // The NodeList we're going apply to
                 const double phi0,                               // Initial porosity
                 const double Pe,                                 // Elastic pressure threshold
                 const double Pt,                                 // Transition pressure (Pe <= Pt)
                 const double Ps,                                 // Solid transition pressure (from fit, Pt <= Ps)
                 const double alphae,                             // distension at Pe
                 const double alphat,                             // Transition distension
                 const double n1,                                 // Fitted exponent for plastic distention evolution
                 const double n2,                                 // Fitted exponent for plastic distention evolution
                 const double cS0,                                // Reference sound speed at full density
                 const double c0,                                 // Reference sound speed at initial porosity
                 const double rhoS0,                              // Reference solid density
                 const bool jutziStateUpdate);                    // Apply state update rules from Jutzi 2008

  PalphaPorosity(const SolidNodeList<Dimension>& nodeList,        // The NodeList we're going apply to
                 const Field<Dimension, Scalar>& phi0,            // Initial porosity
                 const double Pe,                                 // Elastic pressure threshold
                 const double Pt,                                 // Transition pressure (Pe <= Pt)
                 const double Ps,                                 // Solid transition pressure (from fit, Pt <= Ps)
                 const double alphae,                             // distension at Pe
                 const double alphat,                             // Transition distension
                 const double n1,                                 // Fitted exponent for plastic distention evolution
                 const double n2,                                 // Fitted exponent for plastic distention evolution
                 const double cS0,                                // Reference sound speed at full density
                 const Field<Dimension, Scalar>& c0,              // Reference sound speed at initial porosity
                 const double rhoS0,                              // Reference solid density
                 const bool jutziStateUpdate);                    // Apply state update rules from Jutzi 2008

  virtual ~PalphaPorosity();

  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;
  //............................................................................

  //............................................................................
  // Methods required for restarting.
  virtual std::string label() const override { return "PalphaPorosity " + mNodeList.name(); }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //............................................................................

  // Access the material parameters.
  double Pe() const;
  double Pt() const;
  double Ps() const;
  double alphae() const;
  double alphat() const;
  double n1() const;
  double n2() const;
  const Field<Dimension, Scalar>& partialPpartialEps() const;
  const Field<Dimension, Scalar>& partialPpartialRho() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mPe, mPt, mPs, mAlphae, mAlphat, mn1, mn2;
  Field<Dimension, Scalar> mdPdU, mdPdR;

  using PorosityModel<Dimension>::mJutziStateUpdate;
  using PorosityModel<Dimension>::mRhoS0;
  using PorosityModel<Dimension>::mcS0;
  using PorosityModel<Dimension>::mKS0;
  using PorosityModel<Dimension>::mc0;
  using PorosityModel<Dimension>::mMaxAbsDalphaDt;
  using PorosityModel<Dimension>::mNodeList;
  using PorosityModel<Dimension>::mAlpha;
  using PorosityModel<Dimension>::mAlpha0;
  using PorosityModel<Dimension>::mfDS;

  // Disallow default constructor
  PalphaPorosity();
};

}

#include "PalphaPorosityInline.hh"

#endif
