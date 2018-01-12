//------------------------------------------------------------------------------
// Iterate the ideal H algorithm to converge on a new H field.
// This routine replaces the H field in place.
//------------------------------------------------------------------------------
#ifndef __Spheral_iterateIdealH__
#define __Spheral_iterateIdealH__

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "DataBase/DataBase.hh"
#include "Boundary/Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/SmoothingScaleBase.hh"

namespace Spheral {
template<typename Dimension>
void
iterateIdealH(DataBaseSpace::DataBase<Dimension>& dataBase,
              const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
              const KernelSpace::TableKernel<Dimension>& W,
              const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
              const int maxIterations = 100,
              const double tolerance = 1.0e-10,
              const double nPerhForIteration = 0.0,
              const bool sphericalStart = false,
              const bool fixDeterminant = false
              );
}

#endif
