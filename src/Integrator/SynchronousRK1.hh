//---------------------------------Spheral++----------------------------------//
// SynchronousRK1 -- Advance the set of Physics packages in time using first
// order Runge-Kutta -- i.e., forward Euler.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Tue Aug  1 14:44:27 PDT 2006
//----------------------------------------------------------------------------//
#ifndef SynchronousRK1_HH
#define SynchronousRK1_HH

#include "Integrator.hh"

namespace Spheral {

template<typename Dimension>
class SynchronousRK1: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  SynchronousRK1();
  SynchronousRK1(DataBase<Dimension>& dataBase);
  SynchronousRK1(DataBase<Dimension>& dataBase,
                 const std::vector<Physics<Dimension>*>& physicsPackages);

  // Destructor.
  ~SynchronousRK1();

  // Assignment.
  SynchronousRK1& operator=(const SynchronousRK1& rhs);

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const override { return "SynchronousRK1"; }

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
