//---------------------------------Spheral++----------------------------------//
// An experimental hour glass control algorithm for SPH, based on the ASPH
// second moment ideas.
//
// Created by JMO, Sun Jan 15 21:19:53 PST 2006
//----------------------------------------------------------------------------//
#ifndef __Spheral__SecondMomentHourGlassControl__
#define __Spheral__SecondMomentHourGlassControl__

#include "Physics/Physics.hh"

namespace Spheral {

namespace FieldSpace {
  template<typename Dimension, typename DataType> class FieldList;
}
namespace KernelSpace {
  template<typename Dimension> class TableKernel;
}

namespace PhysicsSpace {

template<typename Dimension>
class SecondMomentHourglassControl : 
    public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  SecondMomentHourglassControl(const KernelSpace::TableKernel<Dimension>& W,
                               const double multiplier = 0.05,
                               const double maxAccelerationFactor = 0.001);
                               

  // Destructor.
  virtual ~SecondMomentHourglassControl();

  //******************************************************************************//
  // Methods all Physics packages must provide.
  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBaseSpace::DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBaseSpace::DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
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

  // Local copy of the last acceleration due to this algorithm.
  const FieldSpace::FieldList<Dimension, Vector>& acceleration() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const KernelSpace::TableKernel<Dimension>& mW;
  double mMultiplier;
  double mMaxAccelerationFactor;
#ifndef __GCCXML__
  mutable FieldSpace::FieldList<Dimension, Vector> mAcceleration;
#endif
};

}
}

#ifndef __GCCXML__
#include "SecondMomentHourglassControlInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace PhysicsSpace {
    template<typename Dimension> class SecondMomentHourglassControl;
  }
}

#endif
