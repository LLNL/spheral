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
  using FacetedVolume = typename Dimension::FacetedVolume;

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

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Similarly packages might want a hook to do some post-step finalizations.
  // Really we should rename this post-step finalize.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  // Access our internal data
  const TableKernel<Dimension>&                          WT()            const { return mWT; }
  const FieldList<Dimension, Scalar>&                    zerothMoment()  const { return mZerothMoment; }
  const FieldList<Dimension, Vector>&                    firstMoment()   const { return mFirstMoment; }
  const FieldList<Dimension, SymTensor>&                 secondMoment()  const { return mSecondMoment; }
  const FieldList<Dimension, FacetedVolume>&             cells()         const { return mCells; }
  const FieldList<Dimension, Vector>&                    deltaCentroid() const { return mDeltaCentroid; }
  const FieldList<Dimension, SymTensor>&                 cellSecondMoment() const { return mCellSecondMoment; }

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
  FieldList<Dimension, SymTensor> mSecondMoment, mCellSecondMoment;

  // Voronoi stuff
  FieldList<Dimension, FacetedVolume> mCells;
  FieldList<Dimension, Vector> mDeltaCentroid;
};

}

#endif
