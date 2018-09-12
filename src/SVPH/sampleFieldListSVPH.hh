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

namespace Spheral {

  template<typename Dimension, typename DataType>
  FieldList<Dimension, DataType>
  sampleFieldListSVPH(const FieldList<Dimension, DataType>& fieldList,
                      const FieldList<Dimension, typename Dimension::Vector>& position,
                      const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                      const ConnectivityMap<Dimension>& connectivityMap,
                      const TableKernel<Dimension>& W,
                      const Mesh<Dimension>& mesh,
                      const bool firstOrderConsistent);

}

#endif
