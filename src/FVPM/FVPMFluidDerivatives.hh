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
namespace NodeSpace {
  template<typename Dimension> class FluidNodeList;
}
namespace KernelSpace {
  template<typename Dimension> class TableKernel;
}
namespace FileIOSpace {
  class FileIO;
}

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

namespace FVPMSpace {

template<typename Dimension>
class FVPMFluidDerivatives:
    public NodeSpace::FluidDerivativeProducer<Dimension> {

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
                    const KernelSpace::TableKernel<Dimension>& W,
                    const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                    FieldSpace::Field<Dimension, Scalar>& massDensity) const;

  // The required method for updating the weights.
  virtual
  void
  updateWeight(const State<Dimension>& state,
               const KernelSpace::TableKernel<Dimension>& W,
               const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
               FieldSpace::Field<Dimension, Scalar>& weight) const;

  // Update the correction fields. Unclear what this is in the context of FVPM.
  virtual
  void
  updateCorrections(const State<Dimension>& state,
                    const KernelSpace::TableKernel<Dimension>& W,
                    const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                    FieldSpace::Field<Dimension, Scalar>& omegaGradh,
                    FieldSpace::Field<Dimension, Scalar>& A,
                    FieldSpace::Field<Dimension, Vector>& B,
                    FieldSpace::Field<Dimension, Vector>& C,
                    FieldSpace::Field<Dimension, Tensor>& D,
                    FieldSpace::Field<Dimension, Vector>& gradA,
                    FieldSpace::Field<Dimension, Tensor>& gradB,
                    const bool useGradhCorrections) const;

  // The required method to actually compute the derivatives.
  virtual
  void
  calculateDerivatives(const Scalar time,
                       const Scalar dt,
                       const NodeSpace::FluidNodeList<Dimension>& nodeList,
                       const Scalar nPerh,
                       const bool XSPH,
                       const bool compatibleEnergyEvolution,
                       const Scalar epsTensile,
                       const KernelSpace::TableKernel<Dimension>& W,
                       const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                       const State<Dimension>& state,
                       StateDerivatives<Dimension>& derivatives) const;

protected:
  //--------------------------- Protected Interface ---------------------------//

private:
  //--------------------------- Private Interface ---------------------------//
};

}
}

#endif
