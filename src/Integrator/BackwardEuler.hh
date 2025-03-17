//---------------------------------Spheral++----------------------------------//
// BackwardEuler -- Advance the set of Physics packages in time using first
// order Backward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Mon Oct 21 14:32:05 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_BackwardEuler__
#define __Spheral_BackwardEuler__

#include "Integrator/ImplicitIntegrator.hh"

namespace Spheral {

template<typename Dimension>
class BackwardEuler: public ImplicitIntegrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  BackwardEuler(DataBase<Dimension>& dataBase,
                const std::vector<Physics<Dimension>*> physicsPackages = std::vector<Physics<Dimension>*>(),
                const Scalar beta = 1.0,
                const Scalar ftol = 1.0e-8,
                const Scalar steptol = 1.0e-8,
                const size_t maxIterations = 100u);
  BackwardEuler& operator=(const BackwardEuler& rhs) = default;
  virtual ~BackwardEuler() = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // Access internal state
  Scalar beta()                const { return mBeta; }
  size_t maxIterations()       const { return mMaxIters; }
  Scalar ftol()                const { return mftol; }
  Scalar steptol()             const { return msteptol; }

  void beta(const Scalar x)          { mBeta = x; }
  void maxIterations(const size_t x) { mMaxIters = x; }
  void ftol(const Scalar x)          { mftol = x; }
  void steptol(const Scalar x)       { msteptol = x; }

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const override { return "BackwardEuler"; }

  //--------------------------- Public Interface ---------------------------//
private:
  Scalar mBeta, mftol, msteptol;
  size_t mMaxIters;
};

}

#endif
