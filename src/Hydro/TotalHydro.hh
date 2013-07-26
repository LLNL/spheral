//---------------------------------Spheral++----------------------------------//
// TotalHydro -- The Hydro class for methods that evolve total conserved
// quantities (mass, momentum, and energy) rather than specific quantities.
//
// Created by JMO, Mon Jul 12 21:07:52 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_TotalHydro_hh__
#define __Spheral_TotalHydro_hh__

#include <float.h>
#include <string>
#ifndef __GCCXML__
#include <vector>
#include "DataOutput/registerWithRestart.hh"
#else
#include "fakestl.hh"
#endif

#include "Physics/GenericHydro.hh"
#include "Hydro.hh"  // Needed for enum definitions.

namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  namespace ArtificialViscositySpace {
    template<typename Dimension> class ArtificialViscosity;
  }
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace PhysicsSpace {

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
  TotalHydro(const KernelSpace::TableKernel<Dimension>& W,
             const KernelSpace::TableKernel<Dimension>& WPi,
             ArtificialViscositySpace::ArtificialViscosity<Dimension>& Q,
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
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                     State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs);

  // Post-state update jobs.
  virtual 
  void postStateUpdate(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                       State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivatives) const;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs);

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
  const FieldSpace::FieldList<Dimension, SymTensor>& Hideal() const;
  const FieldSpace::FieldList<Dimension, int>& timeStepMask() const;
  const FieldSpace::FieldList<Dimension, Scalar>& pressure() const;
  const FieldSpace::FieldList<Dimension, Scalar>& soundSpeed() const;
  const FieldSpace::FieldList<Dimension, Scalar>& positionWeight() const;
  const FieldSpace::FieldList<Dimension, Scalar>& weightedNeighborSum() const;
  const FieldSpace::FieldList<Dimension, Scalar>& volume() const;
  const FieldSpace::FieldList<Dimension, Scalar>& totalEnergy() const;
  const FieldSpace::FieldList<Dimension, Vector>& linearMomentum() const;

  const FieldSpace::FieldList<Dimension, Scalar>& DVDt() const;
  const FieldSpace::FieldList<Dimension, Scalar>& DEDt() const;
  const FieldSpace::FieldList<Dimension, Vector>& DpmomDt() const;
  const FieldSpace::FieldList<Dimension, SymTensor>& massSecondMoment() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "Hydro"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  HEvolutionType mHEvolution;
  Scalar mhmin, mhmax, mhratiomin;

  // Some internal scratch fields.
  FieldSpace::FieldList<Dimension, SymTensor> mHideal;
  mutable FieldSpace::FieldList<Dimension, int> mTimeStepMask;
  mutable FieldSpace::FieldList<Dimension, Scalar> mPressure, mSoundSpeed, mPositionWeight, mWeightedNeighborSum, mVolume, mTotalEnergy, mDVDt, mDEDt;
  mutable FieldSpace::FieldList<Dimension, Vector> mLinearMomentum, mDpmomDt;
  mutable FieldSpace::FieldList<Dimension, SymTensor> mMassSecondMoment;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif

  // No default constructor, copying, or assignment.
  TotalHydro();
  TotalHydro(const TotalHydro&);
  TotalHydro& operator=(const TotalHydro&);
};

}
}

#ifndef __GCCXML__
#include "TotalHydroInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class TotalHydro;
  }
}

#endif
