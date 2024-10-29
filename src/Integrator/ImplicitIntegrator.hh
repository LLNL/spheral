//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrator -- Abstract base class for implicit in time integrators.
//
// Created by JMO, Fri Oct 18 16:34:11 PDT 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_ImplicitIntegrator__
#define __Spheral_ImplicitIntegrator__

#include "Integrator/Integrator.hh"

namespace Spheral {

template<typename Dimension>
class ImplicitIntegrator: public Integrator<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PackageIterator = typename std::vector<Physics<Dimension>*>::iterator;
  using ConstPackageIterator = typename std::vector<Physics<Dimension>*>::const_iterator;

  using BoundaryIterator = typename std::vector<Boundary<Dimension>*>::iterator;
  using ConstBoundaryIterator = typename std::vector<Boundary<Dimension>*>::const_iterator;

  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ResidualType = typename Physics<Dimension>::ResidualType;

  // Constructors.
  ImplicitIntegrator(DataBase<Dimension>& dataBase,
                     const Scalar tol = 1.0e-6);
  ImplicitIntegrator(DataBase<Dimension>& dataBase,
                     const std::vector<Physics<Dimension>*>& physicsPackages,
                     const Scalar tol = 1.0e-6);

  // Destructor.
  virtual ~ImplicitIntegrator();

  // Assignment.
  ImplicitIntegrator& operator=(const ImplicitIntegrator& rhs);

  // Find the maximum residual difference in the states
  virtual Scalar computeResiduals(const State<Dimension>& state1,
                                  const State<Dimension>& state0) const;

  // Internal data
  Scalar convergenceTolerance()                  const { return mTol; }
  void convergenceTolerance(const Scalar x)            { mTol = x; }

  // Forbidden methods
  ImplicitIntegrator() = delete;

protected:
  //-------------------------- Protected Interface --------------------------//
  // Override the package dt method to call the implicit version
  virtual TimeStepType dt(const Physics<Dimension>* pkg,
                          const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

private:
  //-------------------------- Private Interface --------------------------//
  Scalar mTol;
};

}

#endif
