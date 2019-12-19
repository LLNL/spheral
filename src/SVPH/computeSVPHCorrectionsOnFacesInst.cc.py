text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/computeSVPHCorrectionsOnFaces.cc"
#include "Geometry/Dimension.hh"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

template 
void
computeSVPHCorrectionsOnFaces(const Mesh<Dim< %(ndim)s > >& mesh,
                              const TableKernel<Dim< %(ndim)s > >& W,
                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& volume,
                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                              const Physics<Dim< %(ndim)s > >::ConstBoundaryIterator& boundaryBegin,
                              const Physics<Dim< %(ndim)s > >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim< %(ndim)s >::Scalar>& A,
                              std::vector<Dim< %(ndim)s >::Vector>& B);

}
"""
