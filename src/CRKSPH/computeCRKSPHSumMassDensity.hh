//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHSumMassDensity__
#define __Spheral__computeCRKSPHSumMassDensity__

#include <vector>

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPHCorrectionParams.hh"

namespace Spheral {
  namespace CRKSPHSpace {

    template<typename Dimension>
    void
    computeCRKSPHSumMassDensity(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                const KernelSpace::TableKernel<Dimension>& W,
                                const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                                const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                                const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& vol,
                                const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                                const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& A,
                                const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& B,
                                const FieldSpace::FieldList<Dimension, typename Dimension::Tensor>& C,
                                const CRKOrder correctionOrder,
                                FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& massDensity);
  }
}

#endif
