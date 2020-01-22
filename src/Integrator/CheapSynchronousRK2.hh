//---------------------------------Spheral++----------------------------------//
// CheapSynchronousRK2 -- Advance the set of Physics packages in time using 
// my standard cheat on the second order Runge-Kutta.  We just skip the end of
// step derivative evaluation, and assume that we can use the previous cycles
// mid step derivatives to determine this steps mid step state.
//
// Created by JMO, Mon Jun 12 17:52:27 PDT 2000
//----------------------------------------------------------------------------//
#ifndef CheapSynchronousRK2_HH
#define CheapSynchronousRK2_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
class CheapSynchronousRK2: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  CheapSynchronousRK2();
  CheapSynchronousRK2(DataBase<Dimension>& dataBase);
  CheapSynchronousRK2(DataBase<Dimension>& dataBase,
                      const std::vector<Physics<Dimension>*>& physicsPackages);

  // Destructor.
  ~CheapSynchronousRK2();

  // Assignment.
  CheapSynchronousRK2& operator=(const CheapSynchronousRK2& rhs);

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const override { return "CheapSynchronousRK2"; }

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class CheapSynchronousRK2;
}

#endif
