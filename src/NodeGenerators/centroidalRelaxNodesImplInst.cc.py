text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "centroidalRelaxNodesImpl.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template
unsigned
centroidalRelaxNodesImpl(DataBaseSpace::DataBase< Dim<%(ndim)s > >& db,
                         const std::vector<Dim< %(ndim)s >::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<Dim< %(ndim)s >::FacetedVolume> >& holes,
                         const KernelSpace::TableKernel<Dim< %(ndim)s > >& W,
                         const PythonBoundFunctors::SpheralFunctor<Dim< %(ndim)s >::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<Dim< %(ndim)s >::Vector, Dim< %(ndim)s >::Vector>& gradrhofunc,
                         const bool useGradRhoFunc,
                         std::vector<BoundarySpace::Boundary<Dim< %(ndim)s > >*>& boundaries,
                         const unsigned maxIterations,
                         const double fracTol,
                         const CRKSPHSpace::CRKOrder correctionOrder,
                         const double centroidFrac,
                         FieldSpace::FieldList<Dim< %(ndim)s >, double>& vol,
                         FieldSpace::FieldList<Dim< %(ndim)s >, int>& surfacePoint,
                         FieldSpace::FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::FacetedVolume>& cells);
}
"""
