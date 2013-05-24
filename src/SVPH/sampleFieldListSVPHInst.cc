//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "sampleFieldListSVPH.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace SVPHSpace {

using std::vector;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using MeshSpace::Mesh;
using BoundarySpace::Boundary;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
// Scalar
template 
FieldList<Dim<1>, Dim<1>::Scalar> 
sampleFieldListSVPH<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& positions,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim<1> >& connectivityMap,
                                            const TableKernel< Dim<1> >& W,
                                            const Mesh<Dim<1> >& mesh,
                                            const vector<Boundary<Dim<1> >*>& boundaries,
                                            const bool firstOrderConsistent);

// Vector
template 
FieldList<Dim<1>, Dim<1>::Vector> 
sampleFieldListSVPH<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& positions,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim<1> >& connectivityMap,
                                            const TableKernel< Dim<1> >& W,
                                            const Mesh<Dim<1> >& mesh,
                                            const vector<Boundary<Dim<1> >*>& boundaries,
                                            const bool firstOrderConsistent);

}
}
