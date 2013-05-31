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

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
// Scalar
template 
FieldList<Dim<1>, Dim<1>::Vector> 
gradientFieldListSVPH<Dim<1>, Dim<1>::Scalar>(const FieldList<Dim<1>, Dim<1>::Scalar>& fieldList,
                                              const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                              const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                              const ConnectivityMap<Dim<1> >& connectivityMap,
                                              const TableKernel< Dim<1> >& W,
                                              const Mesh<Dim<1> >& mesh,
                                              const bool firstOrderConsistent);

// Vector
template 
FieldList<Dim<1>, Dim<1>::Tensor> 
gradientFieldListSVPH<Dim<1>, Dim<1>::Vector>(const FieldList<Dim<1>, Dim<1>::Vector>& fieldList,
                                            const FieldList<Dim<1>, Dim<1>::Vector>& position,
                                            const FieldList<Dim<1>, Dim<1>::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim<1> >& connectivityMap,
                                            const TableKernel< Dim<1> >& W,
                                            const Mesh<Dim<1> >& mesh,
                                            const bool firstOrderConsistent);

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
// Scalar
template 
FieldList<Dim<2>, Dim<2>::Vector> 
gradientFieldListSVPH<Dim<2>, Dim<2>::Scalar>(const FieldList<Dim<2>, Dim<2>::Scalar>& fieldList,
                                              const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                              const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                              const ConnectivityMap<Dim<2> >& connectivityMap,
                                              const TableKernel< Dim<2> >& W,
                                              const Mesh<Dim<2> >& mesh,
                                              const bool firstOrderConsistent);

// Vector
template 
FieldList<Dim<2>, Dim<2>::Tensor> 
gradientFieldListSVPH<Dim<2>, Dim<2>::Vector>(const FieldList<Dim<2>, Dim<2>::Vector>& fieldList,
                                            const FieldList<Dim<2>, Dim<2>::Vector>& position,
                                            const FieldList<Dim<2>, Dim<2>::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim<2> >& connectivityMap,
                                            const TableKernel< Dim<2> >& W,
                                            const Mesh<Dim<2> >& mesh,
                                            const bool firstOrderConsistent);

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
// Scalar
template 
FieldList<Dim<3>, Dim<3>::Vector> 
gradientFieldListSVPH<Dim<3>, Dim<3>::Scalar>(const FieldList<Dim<3>, Dim<3>::Scalar>& fieldList,
                                              const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                              const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                              const ConnectivityMap<Dim<3> >& connectivityMap,
                                              const TableKernel< Dim<3> >& W,
                                              const Mesh<Dim<3> >& mesh,
                                              const bool firstOrderConsistent);

// Vector
template 
FieldList<Dim<3>, Dim<3>::Tensor> 
gradientFieldListSVPH<Dim<3>, Dim<3>::Vector>(const FieldList<Dim<3>, Dim<3>::Vector>& fieldList,
                                            const FieldList<Dim<3>, Dim<3>::Vector>& position,
                                            const FieldList<Dim<3>, Dim<3>::SymTensor>& Hfield,
                                            const ConnectivityMap<Dim<3> >& connectivityMap,
                                            const TableKernel< Dim<3> >& W,
                                            const Mesh<Dim<3> >& mesh,
                                            const bool firstOrderConsistent);

}
}
