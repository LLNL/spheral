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
  PointPotential(double G, double mass, double coreRadius, const Vector& origin);

  // Destructor.
  virtual ~PointPotential();

  // This is the derivative method that all BodyPotential classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  // Provide the timestep appropriate for this package.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  //! Required label for Physics interface.
  virtual std::string label() const { return "PointPotential"; }

  // Get the cumulative potential energy calculated in the last 
  // evaluateDerivatives.
  virtual Scalar extraEnergy() const;

  // The specific potential.
  Scalar specificPotential(const Vector& r) const;

  // Access G.
  Scalar G() const;
  void setG(const Scalar G);

  // Access the mass.
  Scalar mass() const;
  void setMass(const Scalar m);

  // Access the core softening radius.
  Scalar coreRadius() const;
  void setCoreRadius(const Scalar rc);

  // Access the origin.
  const Vector& origin() const;
  void setOrigin(const Vector& origin);

  // The maximum allowed fractional change in a particles potential (for 
  // setting the timestep).
  Scalar deltaPotentialFraction() const;
  void setDeltaPotentialFraction(const Scalar deltaPhi);

private:
  //--------------------------- Public Interface ---------------------------//
  Scalar mG;
  Scalar mMass;
  Scalar mCoreRadius2;
  Vector mOrigin;
  Scalar mDeltaPhiFraction;
  mutable Scalar mPotentialEnergy;

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
