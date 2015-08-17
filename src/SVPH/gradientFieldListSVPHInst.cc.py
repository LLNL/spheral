text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "gradientFieldListSVPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace SVPHSpace {

using std::vector;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using MeshSpace::Mesh;

// Scalar
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector> 
gradientFieldListSVPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& fieldList,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                              const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                              const TableKernel< Dim< %(ndim)s > >& W,
                                              const Mesh<Dim< %(ndim)s > >& mesh,
                                              const bool firstOrderConsistent);

// Vector
template 
FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Tensor> 
gradientFieldListSVPH<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& fieldList,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                                            const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim< %(ndim)s > >& connectivityMap,
                                            const TableKernel< Dim< %(ndim)s > >& W,
                                            const Mesh<Dim< %(ndim)s > >& mesh,
                                            const bool firstOrderConsistent);


}
}
"""
