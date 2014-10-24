//---------------------------------Spheral++----------------------------------//
// Flag the nodes commonly connected within the given H radius as distinct
// fragments.  Result returned as an integer Field of fragement IDs.
//------------------------------------------------------------------------------
#ifndef __Spheral_identifyFragment__
#define __Spheral_identifyFragment__

namespace Spheral {

  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }

  template<typename Dimension>
  FieldSpace::Field<Dimension, int>
  computeFragmentField(const NodeSpace::NodeList<Dimension>& nodeList,
                       const double linkRadius,
                       const FieldSpace::Field<Dimension, typename Dimension::SymTensor>& damage,
                       const double damageThreshold,
                       const bool assignDustToFragments);

}

#endif
