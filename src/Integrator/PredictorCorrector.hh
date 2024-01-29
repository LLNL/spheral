//---------------------------------Spheral++----------------------------------//
// Second order predictor corrector time integrator.
// Based on the predictor corrector scheme described in Monaghans tensile
// instability paper (Monaghan 2000, JCP 159, 290.)
//
// Created by JMO, Wed Dec  4 21:39:36 PST 2002
//----------------------------------------------------------------------------//
#ifndef PredictorCorrector_HH
#define PredictorCorrector_HH

#include "Integrator.hh"

#include <vector>

namespace Spheral {

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
  PredictorCorrector(DataBase<Dimension>& dataBase);
  PredictorCorrector(DataBase<Dimension>& dataBase,
                     const std::vector<Physics<Dimension>*>& physicsPackages);

  // Destructor.
  ~PredictorCorrector();

  // Assignment.
  PredictorCorrector& operator=(const PredictorCorrector& rhs);

  // All Integrators are required to provide the single cycle method.
  virtual bool step(Scalar maxTime,
                    State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) override;

  // We need to make the simpler form of step visible!
  using Integrator<Dimension>::step;

  // Restart methods.
  virtual std::string label() const override { return "PredictorCorrector"; }

};

}

#endif
