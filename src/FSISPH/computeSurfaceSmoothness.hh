//---------------------------------Spheral++------------------------------------
// Quantifies smoothness by comparing the node's interface normal to the normals
// of neighbors belonging to another nodeList. 
//    1.0 -- perfect alignment w/ no "submerged" nodes
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSurfaceSmoothness__
#define __Spheral__computeSurfaceSmoothness__

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;


template<typename Dimension>
void
computeSurfaceSmoothness(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Vector>& surfaceNormals,
                            FieldList<Dimension, typename Dimension::Scalar>& surfaceFraction,
                            FieldList<Dimension, typename Dimension::Scalar>& surfaceNeighborFraction,
                            FieldList<Dimension, typename Dimension::Scalar>& surfaceSmoothness);

}

 
 #endif