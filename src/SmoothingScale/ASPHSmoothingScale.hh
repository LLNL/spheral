//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScale
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 15:01:13 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_ASPHSmooothingScale__
#define __Spheral_ASPHSmooothingScale__

#include "SmoothingScale/SmoothingScaleBase.hh"

namespace Spheral {

template<typename Dimension>
class ASPHSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  ASPHSmoothingScale(const HEvolutionType HUpdate,
                     const TableKernel<Dimension>& W);
  ASPHSmoothingScale() = delete;
  virtual ~ASPHSmoothingScale() {}

  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Access our internal data
  const TableKernel<Dimension>& WT() const;
  const FieldList<Dimension, Scalar>& zerothMoment() const;
  const FieldList<Dimension, Vector>& firstMoment() const;
  const FieldList<Dimension, SymTensor>& secondMoment() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "ASPHSmoothingScale"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mWT;
  FieldList<Dimension, Scalar> mZerothMoment;
  FieldList<Dimension, Vector> mFirstMoment;
  FieldList<Dimension, SymTensor> mSecondMoment;
};

}

#endif
