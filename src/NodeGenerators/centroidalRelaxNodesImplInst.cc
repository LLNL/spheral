//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeGenerators/centroidalRelaxNodesImpl.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template
unsigned
centroidalRelaxNodesImpl(DataBase< Dim<1 > >& db,
                         const std::vector<Dim<1>::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<Dim<1>::FacetedVolume> >& holes,
                         const TableKernel<Dim<1> >& W,
                         const PythonBoundFunctors::SpheralFunctor<Dim<1>::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<Dim<1>::Vector, Dim<1>::Vector>& gradrhofunc,
                         const bool rhoConst,
                         const bool useGradRhoFunc,
                         std::vector<Boundary<Dim<1> >*>& boundaries,
                         const unsigned maxIterations,
                         const double maxFracTol,
                         const double avgFracTol,
                         const RKOrder correctionOrder,
                         const double centroidFrac,
                         FieldList<Dim<1>, double>& vol,
                         FieldList<Dim<1>, int>& surfacePoint,
                         FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells);
#endif

#if defined(SPHERAL_ENABLE_2D)
template
unsigned
centroidalRelaxNodesImpl(DataBase< Dim<2 > >& db,
                         const std::vector<Dim<2>::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<Dim<2>::FacetedVolume> >& holes,
                         const TableKernel<Dim<2> >& W,
                         const PythonBoundFunctors::SpheralFunctor<Dim<2>::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<Dim<2>::Vector, Dim<2>::Vector>& gradrhofunc,
                         const bool rhoConst,
                         const bool useGradRhoFunc,
                         std::vector<Boundary<Dim<2> >*>& boundaries,
                         const unsigned maxIterations,
                         const double maxFracTol,
                         const double avgFracTol,
                         const RKOrder correctionOrder,
                         const double centroidFrac,
                         FieldList<Dim<2>, double>& vol,
                         FieldList<Dim<2>, int>& surfacePoint,
                         FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells);
#endif

#if defined(SPHERAL_ENABLE_3D)
template
unsigned
centroidalRelaxNodesImpl(DataBase< Dim<3 > >& db,
                         const std::vector<Dim<3>::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<Dim<3>::FacetedVolume> >& holes,
                         const TableKernel<Dim<3> >& W,
                         const PythonBoundFunctors::SpheralFunctor<Dim<3>::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<Dim<3>::Vector, Dim<3>::Vector>& gradrhofunc,
                         const bool rhoConst,
                         const bool useGradRhoFunc,
                         std::vector<Boundary<Dim<3> >*>& boundaries,
                         const unsigned maxIterations,
                         const double maxFracTol,
                         const double avgFracTol,
                         const RKOrder correctionOrder,
                         const double centroidFrac,
                         FieldList<Dim<3>, double>& vol,
                         FieldList<Dim<3>, int>& surfacePoint,
                         FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells);
#endif
}