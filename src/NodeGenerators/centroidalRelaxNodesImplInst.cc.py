text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeGenerators/centroidalRelaxNodesImpl.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template
unsigned
centroidalRelaxNodesImpl(DataBase< Dim<%(ndim)s > >& db,
                         const std::vector<Dim< %(ndim)s >::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<Dim< %(ndim)s >::FacetedVolume> >& holes,
                         const TableKernel<Dim< %(ndim)s > >& W,
                         const PythonBoundFunctors::SpheralFunctor<Dim< %(ndim)s >::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<Dim< %(ndim)s >::Vector, Dim< %(ndim)s >::Vector>& gradrhofunc,
                         const bool rhoConst,
                         const bool useGradRhoFunc,
                         std::vector<Boundary<Dim< %(ndim)s > >*>& boundaries,
                         const unsigned maxIterations,
                         const double fracTol,
                         const CRKOrder correctionOrder,
                         const double centroidFrac,
                         FieldList<Dim< %(ndim)s >, double>& vol,
                         FieldList<Dim< %(ndim)s >, int>& surfacePoint,
                         FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::FacetedVolume>& cells);
}
"""
