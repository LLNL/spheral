//---------------------------------Spheral++----------------------------------//
// An experimental hour glass control algorithm based on an estimate of the
// local third moment of the node distribution.
//
// Created by JMO, Thu Apr  2 09:02:00 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral__ThirdMomentHourGlassControl__
#define __Spheral__ThirdMomentHourGlassControl__

#include "Physics/Physics.hh"

namespace Spheral {

template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class TableKernel;

template<typename Dimension>
class ThirdMomentHourglassControl : 
    public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  ThirdMomentHourglassControl(const DataBase<Dimension>& dataBase,
                              const TableKernel<Dimension>& W,
                              const double multiplier = 0.5,
                              const double maxAccelerationFactor = 0.01);

  // Destructor.
  virtual ~ThirdMomentHourglassControl();

  //******************************************************************************//
  // Methods all Physics packages must provide.
  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs);

  // Label
  virtual std::string label() const { return "SecondMomentHourglassControl"; }

  //******************************************************************************//

  // Parameter controlling the maximum allowed acceleration due to the 
  // hourglass control.
  double maxAccelerationFactor() const;
  void maxAccelerationFactor(const double x);

  // Multiplier for the acceleration.
  double multiplier() const;
  void multiplier(const double x);

  // The third moment field.
  const FieldList<Dimension, ThirdRankTensor>& thirdMoment() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mW;
  double mMultiplier;
  double mMaxAccelerationFactor;
  mutable FieldList<Dimension, ThirdRankTensor> mThirdMoment;
};

}

#include "ThirdMomentHourglassControlInline.hh"

#endif
