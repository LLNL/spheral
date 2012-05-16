//------------------------------------------------------------------------------
// Compute the SPH pressure per point that will produce the input accelerations
// given our standard SPH momentum equation.
// We only cover the 3D case here.
//------------------------------------------------------------------------------
#include "HYPRE.h"

#include "Geometry/Dimension.hh"
#include "computeHydrostaticEquilibriumPressure.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"

namespace Spheral {
namespace SPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using NodeSpace::NodeList;

void
computeSPHHydrostaticEquilibriumPressure(const NeighborSpace::ConnectivityMap<Dim<3> >& connectivityMap,
                                         const KernelSpace::TableKernel<Dim<3> >& W,
                                         const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& position,
                                         const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& mass,
                                         const FieldSpace::FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                                         const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& acceleration,
                                         const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& massDensity,
                                         const double tolerance,
                                         const unsigned maxIterations,
                                         FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& pressure) {

  // TAU timers.
  TAU_PROFILE("computeHydrostaticEquilibriumPressure", "", TAU_USER);

  // Pre-conditions.
  const unsigned numNodeLists = pressure.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(acceleration.size() == numNodeLists);
  REQUIRE(massDensity.size() == numNodeLists);
  REQUIRE(position.numGhostNodes() == 0);

  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  // Zero out the result.
  pressure = 0.0;

  // Figure out the offsets for each NodeList.
  vector<unsigned> offset(1, 0);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    offset.push_back(pressure[nodeListi]->numInternalElements());
  }
  CHECK(offset.size() == numNodeLists + 1);

  // Build the sparse matrix that represents the full pressure gradient operator.
  // We have one of these matrix operators for each dimension, hence the 3 vector.
  const unsigned n = pressure.numInternalNodes();
  vector<Eigen::SparseMatrix<double, Eigen::RowMajor> > M(3, Eigen::SparseMatrix<double, Eigen::RowMajor>(n, n));
  unsigned iglobal = 0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = pressure[nodeList]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const unsigned iglobal = offset[nodeListi] + i;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Get the neighbors for this node (in this NodeList).
      const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);
      for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        for (vector<int>::const_iterator jItr = connectivity.begin();
             jItr != connectivity.end();
             ++jItr) {
          const unsigned j = *jItr;
          const unsigned jglobal = offset[nodeListj] + j;

          // State for node j.
          const Vector& rj = position(nodeListj, j);
          const Scalar& mj = mass(nodeListj, j);
          const SymTensor& Hj = H(nodeListj, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel gradient.
          const Vector rji = rj - ri;
          const Vector etai = Hi*rji;
          const Vector etaj = Hj*rji;
          const Scalar etaMagi = etai.magnitude();
          const Scalar etaMagj = etaj.magnitude();
          const Vector gradWij = 0.5*(etai*W.gradValue(etaMagi, Hdeti) +
                                      etaj*W.gradValue(etaMagj, Hdetj));
          M[0](iglobal, jglobal) = mj*gradWij.x();
          M[1](iglobal, jglobal) = mj*gradWij.y();
          M[2](iglobal, jglobal) = mj*gradWij.z();
        }
      }
    }
  }
}

}
}
