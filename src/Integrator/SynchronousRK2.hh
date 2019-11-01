//---------------------------------Spheral++----------------------------------//
// SynchronousRK2 -- Advance the set of Physics packages in time using second
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun  1 22:26:45 PDT 2000
//----------------------------------------------------------------------------//
#ifndef SynchronousRK2_HH
#define SynchronousRK2_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class SynchronousRK2: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  SynchronousRK2();
  SynchronousRK2(DataBase<Dimension>& dataBase);
  SynchronousRK2(DataBase<Dimension>& dataBase,
                 const std::vector<Physics<Dimension>*>& physicsPackages);

  // Destructor.
  ~SynchronousRK2();

  // Assignment.
  SynchronousRK2& operator=(const SynchronousRK2& rhs);

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const override { return "SynchronousRK2"; }

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class SynchronousRK2;
}

#endif
