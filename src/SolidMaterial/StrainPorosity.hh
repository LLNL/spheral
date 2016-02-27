//---------------------------------Spheral++----------------------------------//
// StrainPorosity.hh
//
// An implementation of strain-alpha porosity model described in two papers:
//
// Wunnemann, Collins, & Melosh, Icarus, 180, 514-527 (2006)
// "A strain-based porosity model for use in hydrocode simulations of impacts
//  and implications for transient crater growth in porous targets"
//
// Collins, G. S., Melosh, H. J., & Wünnemann, K. (2011).
// "Improvements to the ɛ-α porous compaction model for simulating impacts into
// high-porosity solar system objects.
// International Journal of Impact Engineering, 38(6), 434–439.
// http://doi.org/10.1016/j.ijimpeng.2010.10.013
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
// StrainPorosity is the physics module which time evolves the distention 
// parameter (alpha) and gives it to the PorousEquationOfState.
//
// Created by JMO, Fri Jun  1 16:16:26 PDT 2012
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrainPorosity__
#define __Spheral_StrainPorosity__

#include <limits>
#include "PorousEquationOfState.hh"
#include "PorousStrengthModel.hh"
#include "Physics/Physics.hh"
#include "DataOutput/registerWithRestart.hh"

namespace Spheral {
namespace SolidMaterial {

template<typename Dimension>
class StrainPorosity: 
    public PhysicsSpace::Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename PhysicsSpace::Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors, destructors.
  StrainPorosity(PorousEquationOfState<Dimension>& porousEOS,     // Porous EOS we're going to modify
                 PorousStrengthModel<Dimension>& porousStrength,  // Porous strength model we're going to modify
                 const NodeSpace::NodeList<Dimension>& nodeList,  // The NodeList we're going apply to
                 const double phi0,                               // Initial porosity
                 const double epsE,                               // Elastic-plastic transition strain
                 const double epsX,                               // Threshold strain between compaction regimes
                 const double kappa,                              // Compaction rate
                 const double gammaS0,                            // Reference gamma at full density
                 const double cS0,                                // Reference sound speed at full density
                 const double c0);                                // Reference sound speed at initial porosity
  virtual ~StrainPorosity();

  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Register our state.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);

  // Do any required one-time initializations on problem start up.
  virtual void initializeProblemStartup(DataBaseSpace::DataBase<Dimension>& dataBase);
  //............................................................................

  //............................................................................
  // Methods required for restarting.
  virtual std::string label() const { return "StrainPorosity " + mNodeList.name(); }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //............................................................................

  // Access the material parameters.
  double phi0() const;
  double alpha0() const;
  double epsE() const;
  double epsX() const;
  double epsC() const;
  double kappa() const;
  double gammaS0() const;
  double cS0() const;
  double c0() const;
  const PorousEquationOfState<Dimension>& porousEOS() const;
  const PorousStrengthModel<Dimension>& porousStrength() const;
  const NodeSpace::NodeList<Dimension>& nodeList() const;
  const FieldSpace::Field<Dimension, Scalar>& alpha() const;
  const FieldSpace::Field<Dimension, Scalar>& DalphaDt() const;
  const FieldSpace::Field<Dimension, Scalar>& strain() const;
  const FieldSpace::Field<Dimension, Scalar>& DstrainDt() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mAlpha0, mEpsE, mEpsX, mKappa, mEpsC, mGammaS0, mcS0, mc0;
  PorousEquationOfState<Dimension>& mPorousEOS;
  PorousStrengthModel<Dimension>& mPorousStrength;
  const NodeSpace::NodeList<Dimension>& mNodeList;
  FieldSpace::Field<Dimension, Scalar> mAlpha,  mDalphaDt, mStrain, mDstrainDt;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;

  // Disallow default constructor
  StrainPorosity();
};

}
}

#include "StrainPorosityInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace SolidMaterial {
    template<typename Dimension> class StrainPorosity;
  }
}

#endif
