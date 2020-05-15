//------------------------------------------------------------------------------
// Implement Lloyd's algorithm for centroidal relaxation of fluid points.
//------------------------------------------------------------------------------
#ifndef __Spheral_centroidalRelaxNodes__
#define __Spheral_centroidalRelaxNodes__

#include "Utilities/Functors.hh"
#include "DataBase/DataBase.hh"
#include "NodeList/FluidNodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "RK/RKCorrectionParams.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
unsigned
centroidalRelaxNodesImpl(DataBase<Dimension>& db,
                         const std::vector<typename Dimension::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<typename Dimension::FacetedVolume> >& holes,
                         const TableKernel<Dimension>& W,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, typename Dimension::Vector>& gradrhofunc,
                         const bool rhoConst,
                         const bool useGradRhoFunc,
                         std::vector<Boundary<Dimension>*>& boundaries,
                         const unsigned maxIterations,
                         const double maxFracTol,
                         const double avgFracTol,
                         const RKOrder correctionOrder,
                         const double centroidFrac,
                         FieldList<Dimension, double>& vol,
                         FieldList<Dimension, int>& surfacePoint,
                         FieldList<Dimension, typename Dimension::FacetedVolume>& cells);

}

#endif
