//---------------------------------Spheral++----------------------------------//
// TotalHydro -- The Hydro class for methods that evolve total conserved
// quantities (mass, momentum, and energy) rather than specific quantities.
//
// Created by JMO, Mon Jul 12 21:07:52 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_TotalHydro_hh__
#define __Spheral_TotalHydro_hh__

#include "DataOutput/registerWithRestart.hh"
#include "Physics/GenericHydro.hh"
#include "Hydro.hh"  // Needed for enum definitions.

#include <float.h>
#include <string>
#include <vector>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class TotalHydro: public GenericHydro<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  TotalHydro(const TableKernel<Dimension>& W,
             const TableKernel<Dimension>& WPi,
             ArtificialViscosity<Dimension>& Q,
             const HEvolutionType HUpdate = IdealH,
             const double hmin = DBL_MIN,
             const double hmax = DBL_MAX,
             const double hratiomin = 0.1);

  // Destructor.
  virtual ~TotalHydro();

  // Hydro's main job is to provide the derivatives for our principle hydrodynamic
  // variables (mass density, velocity, and specific thermal energy).
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // Post-state update jobs.
  virtual 
  void postStateUpdate(const Scalar time, 
                       const Scalar dt,
                       const DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  // Flag to select how we want to evolve the H tensor.
  // the continuity equation.
  HEvolutionType HEvolution() const;
  void HEvolution(const HEvolutionType type);

  // Access the maximum allowed smoothing scale.
  Scalar hmin() const;
  Scalar hmax() const;
  void hmin(const Scalar val);
  void hmax(const Scalar val);

  // Access the minimum allowed ratio of smoothing scales in the h tensor.
  Scalar hratiomin() const;
  void hratiomin(const Scalar val);

  // The state field lists we're maintaining.
  const FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldList<Dimension, int>& timeStepMask() const;
  const FieldList<Dimension, Scalar>& pressure() const;
  const FieldList<Dimension, Scalar>& soundSpeed() const;
  const FieldList<Dimension, Scalar>& positionWeight() const;
  const FieldList<Dimension, Scalar>& weightedNeighborSum() const;
  const FieldList<Dimension, Scalar>& volume() const;
  const FieldList<Dimension, Scalar>& totalEnergy() const;
  const FieldList<Dimension, Vector>& linearMomentum() const;

  const FieldList<Dimension, Scalar>& DVDt() const;
  const FieldList<Dimension, Scalar>& DEDt() const;
  const FieldList<Dimension, Vector>& DpmomDt() const;
  const FieldList<Dimension, SymTensor>& massSecondMoment() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "Hydro"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  HEvolutionType mHEvolution;
  Scalar mhmin, mhmax, mhratiomin;

  // Some internal scratch fields.
  FieldList<Dimension, SymTensor> mHideal;
  mutable FieldList<Dimension, int> mTimeStepMask;
  mutable FieldList<Dimension, Scalar> mPressure, mSoundSpeed, mPositionWeight, mWeightedNeighborSum, mVolume, mTotalEnergy, mDVDt, mDEDt;
  mutable FieldList<Dimension, Vector> mLinearMomentum, mDpmomDt;
  mutable FieldList<Dimension, SymTensor> mMassSecondMoment;

  // The restart registration.
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  TotalHydro();
  TotalHydro(const TotalHydro&);
  TotalHydro& operator=(const TotalHydro&);
};

}

#include "TotalHydroInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class TotalHydro;
}

#endif
