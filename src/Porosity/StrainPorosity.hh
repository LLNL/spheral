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
// Created by JMO, Fri Jun  1 16:16:26 PDT 2012
//----------------------------------------------------------------------------//
#ifndef __Spheral_StrainPorosity__
#define __Spheral_StrainPorosity__

#include "Porosity/PorosityModel.hh"
#include "DataOutput/registerWithRestart.hh"

#include <limits>

namespace Spheral {

template<typename Dimension>
class StrainPorosity: 
    public PorosityModel<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors, destructors.
  StrainPorosity(const SolidNodeList<Dimension>& nodeList,        // The NodeList we're going apply to
                 const double phi0,                               // Initial porosity
                 const double epsE,                               // Elastic-plastic transition strain
                 const double epsX,                               // Threshold strain between compaction regimes
                 const double kappa,                              // Compaction rate
                 const double gammaS0,                            // Reference gamma at full density
                 const double cS0,                                // Reference sound speed at full density
                 const double c0,                                 // Reference sound speed at initial porosity
                 const double rhoS0,                              // Reference solid density
                 const bool jutziStateUpdate);                    // Apply state update rules from Jutzi 2008

  StrainPorosity(const SolidNodeList<Dimension>& nodeList,        // The NodeList we're going apply to
                 const Field<Dimension, Scalar>& phi0,            // Initial porosity
                 const double epsE,                               // Elastic-plastic transition strain
                 const double epsX,                               // Threshold strain between compaction regimes
                 const double kappa,                              // Compaction rate
                 const double gammaS0,                            // Reference gamma at full density
                 const double cS0,                                // Reference sound speed at full density
                 const Field<Dimension, Scalar>& c0,              // Reference sound speed at initial porosity
                 const double rhoS0,                              // Reference solid density
                 const bool jutziStateUpdate);                    // Apply state update rules from Jutzi 2008

  virtual ~StrainPorosity();

  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);
  //............................................................................

  //............................................................................
  // Methods required for restarting.
  virtual std::string label() const { return "StrainPorosity " + mNodeList.name(); }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //............................................................................

  // Access the material parameters.
  double epsE() const;
  double epsX() const;
  double kappa() const;
  double gammaS0() const;
  const Field<Dimension, Scalar>& strain() const;
  const Field<Dimension, Scalar>& DstrainDt() const;

private:
  //--------------------------- Private Interface ---------------------------//
  double mEpsE, mEpsX, mKappa, mGammaS0;
  Field<Dimension, Scalar> mStrain, mDstrainDt;

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
  StrainPorosity();
};

}

#include "StrainPorosityInline.hh"

#endif
