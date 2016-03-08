text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "computeSVPHCorrectionsOnFaces.cc"
#include "Geometry/Dimension.hh"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {
namespace SVPHSpace {

template 
void
computeSVPHCorrectionsOnFaces(const MeshSpace::Mesh<Dim< %(ndim)s > >& mesh,
                              const KernelSpace::TableKernel<Dim< %(ndim)s > >& W,
                              const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>& volume,
                              const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position,
                              const FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>& H,
                              const PhysicsSpace::Physics<Dim< %(ndim)s > >::ConstBoundaryIterator& boundaryBegin,
                              const PhysicsSpace::Physics<Dim< %(ndim)s > >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim< %(ndim)s >::Scalar>& A,
                              std::vector<Dim< %(ndim)s >::Vector>& B);

}
}
"""
