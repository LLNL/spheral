//---------------------------------Spheral++----------------------------------//
// Verlet -- Advance the set of Physics packages in time using the second
// order Verlet algorithm.  This method is symplectic in the absence of 
// dissipation.
//
// Based on the description in 
// Monaghan JJ. Smoothed particle hydrodynamics. Reports on progress in physics. 2005
//
// Created by JMO, Sat Aug 23 10:33:33 PDT 2014
//----------------------------------------------------------------------------//
#ifndef Verlet_HH
#define Verlet_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class Verlet: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  Verlet();
  Verlet(DataBase<Dimension>& dataBase);
  Verlet(DataBase<Dimension>& dataBase,
                 const std::vector<Physics<Dimension>*>& physicsPackages);

  // Destructor.
  ~Verlet();

  // Assignment.
  Verlet& operator=(const Verlet& rhs);

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const override { return "Verlet"; }

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class Verlet;
}

#endif
