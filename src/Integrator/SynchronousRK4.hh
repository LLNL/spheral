//---------------------------------Spheral++----------------------------------//
// SynchronousRK4 -- Advance the set of Physics packages in time using fourth
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun 14 10:22:47 PDT 2012
//----------------------------------------------------------------------------//
#ifndef SynchronousRK4_HH
#define SynchronousRK4_HH

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "Integrator.hh"

namespace Spheral {
namespace IntegratorSpace {

template<typename Dimension>
class SynchronousRK4: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  SynchronousRK4();
  SynchronousRK4(DataBaseSpace::DataBase<Dimension>& dataBase);
  SynchronousRK4(DataBaseSpace::DataBase<Dimension>& dataBase,
                 const std::vector<PhysicsSpace::Physics<Dimension>*>& physicsPackages);

  // Destructor.
  ~SynchronousRK4();

  // Assignment.
  SynchronousRK4& operator=(const SynchronousRK4& rhs);

  // All Integrators are required to provide the single cycle method.
  virtual void step(Scalar maxTime);

  // Restart methods.
  virtual std::string label() const { return "SynchronousRK4"; }

private:
  //--------------------------- Private Interface ---------------------------//
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace IntegratorSpace {
    template<typename Dimension> class SynchronousRK4;
  }
}

#endif
