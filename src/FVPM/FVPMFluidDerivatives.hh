//---------------------------------Spheral++----------------------------------//
// FVPMFluidDerivatives
//
// Compute the SPH definition for the fluid derivatives.
//
// Created by JNJ, Sat Jul 10 12:08:50 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_FVPMSpace_FVPMFluidDerivatives__
#define __Spheral_FVPMSpace_FVPMFluidDerivatives__

#include "NodeList/FluidDerivativeProducer.hh"

namespace Spheral {

// Forward declarations.
template<typename Dimension> class FluidNodeList;
template<typename Dimension> class TableKernel;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
class FileIO;

template<typename Dimension>
class FVPMFluidDerivatives:
    public FluidDerivativeProducer<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors, destructor.
  FVPMFluidDerivatives();
  FVPMFluidDerivatives(const FVPMFluidDerivatives& rhs);
  FVPMFluidDerivatives& operator=(const FVPMFluidDerivatives& rhs);
  virtual ~FVPMFluidDerivatives();

  // The required summed definition of the mass density.
  virtual
  void
  updateMassDensity(const State<Dimension>& state,
                    const TableKernel<Dimension>& W,
                    const ConnectivityMap<Dimension>& connectivityMap,
                    Field<Dimension, Scalar>& massDensity) const;

  // The required method for updating the weights.
  virtual
  void
  updateWeight(const State<Dimension>& state,
               const TableKernel<Dimension>& W,
               const ConnectivityMap<Dimension>& connectivityMap,
               Field<Dimension, Scalar>& weight) const;

  // Update the correction fields. Unclear what this is in the context of FVPM.
  virtual
  void
  updateCorrections(const State<Dimension>& state,
                    const TableKernel<Dimension>& W,
                    const ConnectivityMap<Dimension>& connectivityMap,
                    Field<Dimension, Scalar>& omegaGradh,
                    Field<Dimension, Scalar>& A,
                    Field<Dimension, Vector>& B,
                    Field<Dimension, Vector>& C,
                    Field<Dimension, Tensor>& D,
                    Field<Dimension, Vector>& gradA,
                    Field<Dimension, Tensor>& gradB,
                    const bool useGradhCorrections) const;

  // The required method to actually compute the derivatives.
  virtual
  void
  calculateDerivatives(const Scalar time,
                       const Scalar dt,
                       const FluidNodeList<Dimension>& nodeList,
                       const Scalar nPerh,
                       const bool XSPH,
                       const bool compatibleEnergyEvolution,
                       const Scalar epsTensile,
                       const TableKernel<Dimension>& W,
                       const ConnectivityMap<Dimension>& connectivityMap,
                       const State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) const;

protected:
  //--------------------------- Protected Interface ---------------------------//

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#endif
