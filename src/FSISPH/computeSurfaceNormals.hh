//---------------------------------Spheral++------------------------------------
// Compute normals at an interface
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSurfaceNormals__
#define __Spheral__computeSurfaceNormals__

namespace Spheral {

  // Forward declarations.
  template<typename Dimension> class ConnectivityMap;
  template<typename Dimension> class TableKernel;
  template<typename Dimension, typename DataType> class FieldList;


template<typename Dimension>
void
computeSurfaceNormals(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::Scalar>& massDensity,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Vector>& interfaceNormals);

}

 
 #endif