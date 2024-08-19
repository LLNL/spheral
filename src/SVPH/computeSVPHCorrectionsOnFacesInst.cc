//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "SVPH/computeSVPHCorrectionsOnFaces.cc"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

template
void
computeSVPHCorrectionsOnFaces(const Mesh<Dim<1> >& mesh,
                              const TableKernel<Dim<1> >& W,
                              const FieldList<Dim<1>, Dim<1>::Scalar>& volume,
                              const FieldList<Dim<1>, Dim<1>::Vector>& position,
                              const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                              const Physics<Dim<1> >::ConstBoundaryIterator& boundaryBegin,
                              const Physics<Dim<1> >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim<1>::Scalar>& A,
                              std::vector<Dim<1>::Vector>& B);

#endif

#if defined(SPHERAL_ENABLE_2D)

template
void
computeSVPHCorrectionsOnFaces(const Mesh<Dim<2> >& mesh,
                              const TableKernel<Dim<2> >& W,
                              const FieldList<Dim<2>, Dim<2>::Scalar>& volume,
                              const FieldList<Dim<2>, Dim<2>::Vector>& position,
                              const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                              const Physics<Dim<2> >::ConstBoundaryIterator& boundaryBegin,
                              const Physics<Dim<2> >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim<2>::Scalar>& A,
                              std::vector<Dim<2>::Vector>& B);

#endif

#if defined(SPHERAL_ENABLE_3D)

template
void
computeSVPHCorrectionsOnFaces(const Mesh<Dim<3> >& mesh,
                              const TableKernel<Dim<3> >& W,
                              const FieldList<Dim<3>, Dim<3>::Scalar>& volume,
                              const FieldList<Dim<3>, Dim<3>::Vector>& position,
                              const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                              const Physics<Dim<3> >::ConstBoundaryIterator& boundaryBegin,
                              const Physics<Dim<3> >::ConstBoundaryIterator& boundaryEnd,
                              std::vector<Dim<3>::Scalar>& A,
                              std::vector<Dim<3>::Vector>& B);

#endif
}