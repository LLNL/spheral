//---------------------------------Spheral++----------------------------------//
// PolyGravity -- Solve gravity on a polytope (2D or 3D).
//
// Based on Jason Pearl's approximate gravity model.  Currently 2D polygons
// are not really correct, as they use the 3D implementation and therefore
// don't properly reprsent the 2D infinite-rod logarithmic potential.
//
// Created by JMO, Fri Sep 16 15:04:54 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_PolyGravity__
#define __Spheral_PolyGravity__

#include "Gravity/ApproximatePolyhedralGravityModel.hh"
#include "Gravity/TreeGravity.hh"      // Just for GravityTimeStepType enum
#include "Geometry/Dimension.hh"
#include "Physics/GenericBodyForce.hh"
#include "Field/FieldList.hh"

#include <memory>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
class FileIO;

template<typename Dimension>
class PolyGravity: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using Polytope = typename Dimension::FacetedVolume;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  //! Constructor.
  //! \param poly -- polytope we're going to solve on.
  //! \param G -- the gravitational constant.
  //! \param mass -- mass of the polytope
  //! \param ftimestep -- safety factor in [0,1] in setting time steps.
  //! \param timeStepChoice -- algorithm to use choosing time steps.
  PolyGravity(const Polytope& poly,
              const double G,
              const double mass,
              const double ftimestep,
              const GravityTimeStepType timeStepChoice);

  //! Destructor.
  virtual ~PolyGravity();

  //! We augment the generic body force state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  //! This is the derivative method that all BodyForce classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar /*time*/,
                           const Scalar /*dt*/,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  //! Vote on the timestep.  This uses a velocity-limiting rule.
  virtual TimeStepType dt(const DataBase<Dimension>& /*dataBase*/, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& /*derivs*/,
                          const Scalar /*currentTime*/) const;

  //! Initializations on problem start up.
  virtual void initializeProblemStartup(DataBase<Dimension>& db);

  //! This package opts out of building connectivity.
  virtual bool requireConnectivity() const { return false; }

  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const;

  //! Return the gravitational potential created by the particle distribution.
  const FieldList<Dimension, Scalar>& potential() const;

  //! The surface model
  const Dim<3>::FacetedVolume& poly() const;

  //! The gravitational constant we're using.
  double G() const;

  //! The mass of the polytope
  double mass() const;

  //! The current time step scaling factor.
  double ftimestep() const;
  void ftimestep(double x);

  //! The algorithmic choice for setting the time step.
  GravityTimeStepType timeStepChoice() const;
  void timeStepChoice(GravityTimeStepType x);

  //! The dynamical time
  double dynamicalTime() const;

  //! Gravity solver
  const ApproximatePolyhedralGravityModel& solver() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "PolyGravity"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  // Private data.
  double mG, mMass, mftimestep, mDynamicalTime;
  GravityTimeStepType mTimeStepChoice;
  Dim<3>::FacetedVolume mPoly;
  std::shared_ptr<ApproximatePolyhedralGravityModel> mSolver;

  // The potential fields filled in during evaluateDerivates.
  mutable FieldList<Dimension, Scalar> mPotential;
  mutable Scalar mExtraEnergy;

  // Data we need for computing time steps.
  mutable Scalar mDtMinAcc;
  
  // The restart registration.
  RestartRegistrationType mRestart;

  // Default constructor -- disabled.
  PolyGravity();

  // Copy constructor -- disabled.
  PolyGravity(const PolyGravity&);

  // Assignment operator -- disabled.
  PolyGravity& operator=(const PolyGravity&);
};

// Declare explicit specializations.
template<> PolyGravity<Dim<2>>::PolyGravity(const PolyGravity<Dim<2>>::Polytope&, const double, const double, const double, const GravityTimeStepType);
template<> PolyGravity<Dim<3>>::PolyGravity(const PolyGravity<Dim<3>>::Polytope&, const double, const double, const double, const GravityTimeStepType);

}

#endif
