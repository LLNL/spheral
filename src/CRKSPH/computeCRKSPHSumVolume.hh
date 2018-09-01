//---------------------------------Spheral++------------------------------------
// Compute the CRKSPH mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeCRKSPHSumVolume__
#define __Spheral__computeCRKSPHSumVolume__

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
void
computeCRKSPHSumVolume(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                            const KernelSpace::TableKernel<Dimension>& W,
                            const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& vol);

}

#endif
