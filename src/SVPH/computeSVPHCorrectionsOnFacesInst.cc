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
computeSVPHCorrectionsOnFaces(const MeshSpace::Mesh<Dim<1> >& mesh,
                              const KernelSpace::TableKernel<Dim<1> >& W,
                              const FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& volume,
                              const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& position,
                              const FieldSpace::FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                              const PhysicsSpace::Physics<Dim<1> >::ConstBoundaryIterator& boundaryBegin,
                              const PhysicsSpace::Physics<Dim<1> >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim<1>::Scalar>& A,
                              std::vector<Dim<1>::Vector>& B);

template 
void
computeSVPHCorrectionsOnFaces(const MeshSpace::Mesh<Dim<2> >& mesh,
                              const KernelSpace::TableKernel<Dim<2> >& W,
                              const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>& volume,
                              const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& position,
                              const FieldSpace::FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                              const PhysicsSpace::Physics<Dim<2> >::ConstBoundaryIterator& boundaryBegin,
                              const PhysicsSpace::Physics<Dim<2> >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim<2>::Scalar>& A,
                              std::vector<Dim<2>::Vector>& B);

template 
void
computeSVPHCorrectionsOnFaces(const MeshSpace::Mesh<Dim<3> >& mesh,
                              const KernelSpace::TableKernel<Dim<3> >& W,
                              const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& volume,
                              const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& position,
                              const FieldSpace::FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                              const PhysicsSpace::Physics<Dim<3> >::ConstBoundaryIterator& boundaryBegin,
                              const PhysicsSpace::Physics<Dim<3> >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim<3>::Scalar>& A,
                              std::vector<Dim<3>::Vector>& B);

}
}
