//------------------------------------------------------------------------------
// Use SVPH to take the gradient of a FieldList.
//------------------------------------------------------------------------------
#ifndef __Spheral__gradientFieldListSVPH__
#define __Spheral__gradientFieldListSVPH__

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Mesh/Mesh.hh"
#include "Geometry/MathTraits.hh"

#include <vector>

namespace Spheral {

  template<typename Dimension, typename DataType>
  FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
  gradientFieldListSVPH(const FieldList<Dimension, DataType>& fieldList,
                        const FieldList<Dimension, typename Dimension::Vector>& position,
                        const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
                        const ConnectivityMap<Dimension>& connectivityMap,
                        const TableKernel<Dimension>& W,
                        const Mesh<Dimension>& mesh,
                        const bool firstOrderConsistent);

}

#endif
