//---------------------------------Spheral++----------------------------------//
// PointPotential -- Impose a potential from a point mass.
//
// Created by JMO, Sun Mar 30 22:08:55 PST 2003
//----------------------------------------------------------------------------//
#ifndef PointPotential_HH
#define PointPotential_HH

#include "Physics/GenericBodyForce.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class NodeList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class PointPotential: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  PointPotential(double G, double mass, double coreRadius, const Vector origin,
                 const Tensor metric);

  // Destructor.
  virtual ~PointPotential();

  // We augment the generic body force state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // This is the derivative method that all BodyPotential classes must provide.
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivs) const;

  // Provide the timestep appropriate for this package.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  //! Initializations on problem start up.
  virtual void initializeProblemStartup(DataBase<Dimension>& db);

  // Get the cumulative potential energy calculated in the last 
  // evaluateDerivatives.
  virtual Scalar extraEnergy() const;

  // The specific potential.
  Scalar specificPotential(const Vector& r) const;

  // Access G.
  Scalar G() const;
  void G(Scalar G);

  // Access the mass.
  Scalar mass() const;
  void mass(Scalar m);

  // Access the core softening radius.
  Scalar coreRadius() const;
  void coreRadius(Scalar rc);

  // Access the origin.
  const Vector& origin() const;
  void origin(const Vector& origin);

  // The metric for modifying distance calculation
  const Tensor& metric() const;
  void metric(const Tensor& metric);

  //! The current time step scaling factor.
  Scalar ftimestep() const;
  void ftimestep(Scalar x);

  // Return the gravitational potential created by the particle distribution.
  const FieldList<Dimension, Scalar>& potential() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "PointPotential"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Public Interface ---------------------------//
  Scalar mG, mMass, mCoreRadius2, mftimestep;
  Vector mOrigin;
  Tensor mMetric;
  mutable Scalar mDtMinAcc, mTotalPotentialEnergy;
  mutable FieldList<Dimension, Scalar> mPotential;

  // No default constructor, copying, or assignment.
  PointPotential();
  PointPotential(const PointPotential& rhs);
  PointPotential& operator=(const PointPotential& rhs);
};

}

#include "PointPotentialInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class PointPotential;
}

#endif
