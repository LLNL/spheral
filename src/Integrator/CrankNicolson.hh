//---------------------------------Spheral++----------------------------------//
// CrankNicolson -- Advance the set of Physics packages in time using the second
// order Crank-Nicolson method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Wed Mar 19 15:49:24 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_CrankNicolson__
#define __Spheral_CrankNicolson__

#include "Integrator/ImplicitIntegrator.hh"
#include <vector>

namespace Spheral {

template<typename Dimension>
class CrankNicolson: public ImplicitIntegrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  CrankNicolson(DataBase<Dimension>& dataBase,
                const std::vector<Physics<Dimension>*> physicsPackages = std::vector<Physics<Dimension>*>(),
                const Scalar alpha = 0.5,
                const Scalar tol = 1.0e-8,
                const size_t maxIterations = 100u);
  CrankNicolson& operator=(const CrankNicolson& rhs) = default;
  virtual ~CrankNicolson() = default;

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // Access internal state
  Scalar alpha()                          const { return mAlpha; }
  size_t maxIterations()                  const { return mMaxIterations; }
  size_t numExplicitSteps()               const { return mNumExplicitSteps; }
  size_t numImplicitSteps()               const { return mNumImplicitSteps; }

  void alpha(const Scalar x)                    { mAlpha = x; }
  void maxIterations(const size_t x)            { mMaxIterations = x; }

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const override   { return "CrankNicolson"; }

  //--------------------------- Public Interface ---------------------------//
private:
  Scalar mAlpha;
  size_t mMaxIterations, mNumExplicitSteps, mNumImplicitSteps;

  using Integrator<Dimension>::mDtMultiplier;
};

}

#endif
