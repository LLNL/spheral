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
#include "SmoothingScale/ASPHSmoothingScaleUserFilter.hh"
#include "SmoothingScale/ASPHRadialFunctor.hh"
#include "Utilities/Functors.hh"

#include <memory>   // std::shared_ptr
#include <utility>  // std::pair
#include <tuple>

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
  using HidealFilterType = ASPHSmoothingScaleUserFilter<Dimension>;
  using RadialFunctorType = ASPHRadialFunctor<Dimension>;

  // Constructors, destructor.
  ASPHSmoothingScale(const HEvolutionType HUpdate,
                     const TableKernel<Dimension>& W,
                     const bool fixShape = false,
                     const bool radialOnly = false);
  ASPHSmoothingScale() = delete;
  virtual ~ASPHSmoothingScale() {}

  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

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

  // Optional hook to be called at the beginning of a time step.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

 // Similarly packages might want a hook to do some post-step finalizations.
  // Really we should rename this post-step finalize.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  // We require the Voronoi-like cells per point
  virtual bool requireVoronoiCells() const override                               { return (this->HEvolution() == HEvolutionType::IdealH and
                                                                                            not (mFixShape or mRadialOnly)); }

  // Access our internal data
  const TableKernel<Dimension>&                          WT()            const    { return mWT; }
  const FieldList<Dimension, Scalar>&                    zerothMoment()  const    { return mZerothMoment; }
  const FieldList<Dimension, SymTensor>&                 secondMoment()  const    { return mSecondMoment; }
  const FieldList<Dimension, SymTensor>&                 cellSecondMoment() const { return mCellSecondMoment; }
  const FieldList<Dimension, Scalar>&                    radius0() const          { return mRadius0; }

  // Special evolution flags
  bool fixShape() const                                                           { return mFixShape; }
  bool radialOnly() const                                                         { return mRadialOnly; }
  void fixShape(const bool x)                                                     { mFixShape = x; }
  void radialOnly(const bool x)                                                   { mRadialOnly = x; }

  // Optional user functor to manipulate the final ideal H vote
  std::shared_ptr<HidealFilterType> HidealFilter() const                          { return mHidealFilterPtr; }
  void HidealFilter(std::shared_ptr<HidealFilterType> functorPtr)                 { mHidealFilterPtr = functorPtr; }

  // Optional user functor to override the radial unit normal and radius for radialOnly mode
  std::shared_ptr<RadialFunctorType> RadialFunctor() const                        { return mRadialFunctorPtr; }
  void RadialFunctor(std::shared_ptr<RadialFunctorType> functorPtr)               { mRadialFunctorPtr = functorPtr; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "ASPHSmoothingScale"; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mWT;
  FieldList<Dimension, Scalar> mZerothMoment;
  FieldList<Dimension, SymTensor> mSecondMoment, mCellSecondMoment;
  FieldList<Dimension, Scalar> mRadius0;
  std::shared_ptr<HidealFilterType> mHidealFilterPtr;
  std::shared_ptr<RadialFunctorType> mRadialFunctorPtr;
  bool mFixShape, mRadialOnly;
};

}

#endif
