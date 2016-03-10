//---------------------------------Spheral++----------------------------------//
// Second order predictor corrector time integrator.
// Based on the predictor corrector scheme described in Monaghans tensile
// instability paper (Monaghan 2000, JCP 159, 290.)
//
// Created by JMO, Wed Dec  4 21:39:36 PST 2002
//----------------------------------------------------------------------------//
#ifndef PredictorCorrector_HH
#define PredictorCorrector_HH

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "Integrator.hh"

namespace Spheral {
namespace IntegratorSpace {

template<typename Dimension>
class PredictorCorrector: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  PredictorCorrector();
  PredictorCorrector(DataBaseSpace::DataBase<Dimension>& dataBase);
  PredictorCorrector(DataBaseSpace::DataBase<Dimension>& dataBase,
                     const std::vector<PhysicsSpace::Physics<Dimension>*>& physicsPackages);

  // Destructor.
  ~PredictorCorrector();

  // Assignment.
  PredictorCorrector& operator=(const PredictorCorrector& rhs);

  // All Integrators are required to provide the single cycle method.
  virtual void step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs);

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const { return "PredictorCorrector"; }

};
}
}

#else

// Forward declaration.
namespace Spheral {
  namespace IntegratorSpace {
    template<typename Dimension> class PredictorCorrector;
  }
}

#endif
