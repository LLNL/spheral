//---------------------------------Spheral++----------------------------------//
// An experimental hour glass control algorithm for SPH, based on the using
// center of mass estimates in Voronoi cells.
//
// Created by JMO, Tue Jun 28 14:54:03 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral__VoronoiHourGlassControl__
#define __Spheral__VoronoiHourGlassControl__

#include "Physics/Physics.hh"

#include <vector>

namespace Spheral {

template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class TableKernel;

template<typename Dimension>
class VoronoiHourglassControl : 
    public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors.
  // order:    0 => constant in cell
  //           1 => linear slope in cell
  // limiter:  0 => unlimited 
  //           1 => scalar limiter
  //           2 => tensor limiter
  // fraction: \in [0,1], allowed fractional displacement.
  // mask:     Optional mask indicating nodes to filter or not
  //             0 => don't apply filtering.
  //             1 => filter
  VoronoiHourglassControl(const TableKernel<Dimension>& W,
                          const unsigned order,
                          const unsigned limiter,
                          const double fraction,
                          const FieldList<Dimension, int>& mask);

  // Destructor.
  virtual ~VoronoiHourglassControl();

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

  // Finalize method, where we override the positions.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs);

  // Label
  virtual std::string label() const { return "SecondMomentHourglassControl"; }

  //******************************************************************************//
  // The order for the density fit in a cell:
  //   0 => constant
  //   1 => linear
  unsigned order() const;
  void order(const unsigned x);

  // The limiter to be used when finding the density slop in a cell:
  //   0 => unlimited
  //   1 => scalar limiter
  //   2 => tensor limiter
  unsigned limiter() const;
  void limiter(const unsigned x);

  // The allowed fractional displacement.
  double fraction() const;
  void fraction(const double x);

  // Mask indicating which nodes are allowed filtering.
  const FieldList<Dimension, int>& mask() const;
  void mask(const FieldList<Dimension, int>& x);

  // The Kernel.
  const TableKernel<Dimension>& kernel() const;

  // Last gradient of the mass density.
  const FieldList<Dimension, Vector>& gradRho() const;

  // CRKSPH correction fields.
  const FieldList<Dimension, Scalar>& A() const;
  const FieldList<Dimension, Vector>& B() const;
  const FieldList<Dimension, Vector>& C() const;
  const FieldList<Dimension, Tensor>& D() const;
  const FieldList<Dimension, Vector>& gradA() const;
  const FieldList<Dimension, Tensor>& gradB() const;
  const FieldList<Dimension, Scalar>& weight() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mW;
  unsigned mOrder, mLimiter;
  double mFraction;
  FieldList<Dimension, int> mMask;
  FieldList<Dimension, Scalar> mA, mWeight;
  FieldList<Dimension, Vector> mGradRho, mB, mC, mGradA;
  FieldList<Dimension, Tensor> mD, mGradB;
};

}

#include "VoronoiHourglassControlInline.hh"

#endif
