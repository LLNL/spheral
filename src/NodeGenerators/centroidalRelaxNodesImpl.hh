//------------------------------------------------------------------------------
// Implement Lloyd's algorithm for centroidal relaxation of fluid points.
//------------------------------------------------------------------------------
#ifndef __Spheral_centroidalRelaxNodes__
#define __Spheral_centroidalRelaxNodes__

#include <vector>

#include "Utilities/Functors.hh"
#include "DataBase/DataBase.hh"
#include "NodeList/FluidNodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"
#include "CRKSPH/CRKSPHCorrectionParams.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"

namespace Spheral {

template<typename Dimension>
unsigned
centroidalRelaxNodesImpl(DataBaseSpace::DataBase<Dimension>& db,
                         const std::vector<typename Dimension::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<typename Dimension::FacetedVolume> >& holes,
                         const KernelSpace::TableKernel<Dimension>& W,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, typename Dimension::Vector>& gradrhofunc,
                         const bool useGradRhoFunc,
                         std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                         const unsigned maxIterations,
                         const double fracTol,
                         const CRKSPHSpace::CRKOrder correctionOrder,
                         const double centroidFrac,
                         FieldSpace::FieldList<Dimension, double>& vol,
                         FieldSpace::FieldList<Dimension, int>& surfacePoint,
                         FieldSpace::FieldList<Dimension, typename Dimension::FacetedVolume>& cells);

}

#endif
