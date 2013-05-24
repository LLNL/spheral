//------------------------------------------------------------------------------
// Use SVPH to sample a FieldList.
//------------------------------------------------------------------------------
#ifndef __Spheral__sampleFieldListSVPH__
#define __Spheral__sampleFieldListSVPH__

#include <vector>
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Mesh/Mesh.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {
  namespace SVPHSpace {

    template<typename Dimension, typename DataType>
    FieldList<Dimension, DataType>
    sampleFieldListSVPH(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
                        const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& positions,
                        const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                        const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                        const KernelSpace::TableKernel<Dimension>& W,
                        const MeshSpace::Mesh<Dimension>& mesh,
                        const std::vector<BoundarySpace::Boundary<Dimension>*> boundaries,
                        const bool firstOrderConsistent);
  }
}

#endif
