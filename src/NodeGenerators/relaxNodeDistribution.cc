//------------------------------------------------------------------------------
// Centroidally relax a node distribution in a boundary.
// Optionally the user can specify a weighting function for the nodes.
//------------------------------------------------------------------------------
#include <ctime>
#include "relaxNodeDistribution.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"

#ifdef USE_MPI
#include "mpi.h"
#include "Distributed/Communicator.hh"
#endif

namespace Spheral {

using namespace std;

using DataBaseSpace::DataBase;
using BoundarySpace::Boundary;
using KernelSpace::TableKernel;
using FieldSpace::FieldList;
using FieldSpace::Field;
using NeighborSpace::ConnectivityMap;
using NodeSpace::SmoothingScaleBase;

template<typename Dimension>
void
relaxNodeDistribution(DataBaseSpace::DataBase<Dimension>& dataBase,
                      const typename Dimension::FacetedVolume& boundary,
                      const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                      const WeightingFunctor<Dimension>& weighting,
                      const int maxIterations = 100,
                      const double tolerance = 1.0e-10) {
}

}

