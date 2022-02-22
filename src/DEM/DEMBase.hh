//---------------------------------Spheral++----------------------------------//
// DEMBase -- basic DEM package for Spheral++.
//----------------------------------------------------------------------------//
#ifndef __Spheral_DEMBase_hh__
#define __Spheral_DEMBase_hh__

#include <string>
#include "Physics/Physics.hh"
#include "DEM/DEMDimension.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class ContactModelBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class DEMBase: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename DEMDimension<Dimension>::AngularVector RotationType;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  DEMBase(const DataBase<Dimension>& dataBase,
          const Scalar stepsPerCollision,
          const Vector& xmin,
          const Vector& xmax);

  // Destructor.
  virtual ~DEMBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  void initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) override;
                       

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;
  

  // Optionally we can provide a bounding box for use generating the mesh
  // for the Voronoi mass density update.
  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  Scalar stepsPerCollision() const;
  void   stepsPerCollision(Scalar x);

  // access for fieldLists
  const FieldList<Dimension, int>&    timeStepMask() const;
  const FieldList<Dimension, Vector>& DxDt() const;
  const FieldList<Dimension, Vector>& DvDt() const;
  const FieldList<Dimension, RotationType>& omega() const;
  const FieldList<Dimension, RotationType>& DomegaDt() const;

  // inlined and specialized for different dimensions
  virtual Scalar momentOfInertia(const Scalar massi,
                                 const Scalar particleRadiusi) const;

  // // inlined and specialized for different dimensions
  // virtual Vector cross(const Vector contactDirection,
  //                      const RotationType angularVelocity) const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "DEMBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//

  // number of steps per collision time-scale
  Scalar mStepsPerCollision;              

  // Optional bounding box for generating the mesh.
  Vector mxmin, mxmax;

  // Some internal scratch fields.
  FieldList<Dimension, int>      mTimeStepMask;
  FieldList<Dimension, Vector>   mDxDt;
  FieldList<Dimension, Vector>   mDvDt;
  FieldList<Dimension, RotationType> mOmega;
  FieldList<Dimension, RotationType> mDomegaDt;

  // The restart registration.
  RestartRegistrationType mRestart;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  DEMBase();
  DEMBase(const DEMBase&);
  DEMBase& operator=(const DEMBase&);
};

}

#include "DEMBaseInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class DEMBase;
}

#endif
