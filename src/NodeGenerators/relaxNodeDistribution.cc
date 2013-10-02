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
using MeshSpace::Mesh;

template<typename Dimension>
void
relaxNodeDistribution(DataBaseSpace::DataBase<Dimension>& dataBase,
                      const typename Dimension::FacetedVolume& boundary,
                      const std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                      const KernelSpace::TableKernel<Dimension>& W,
                      const NodeSpace::SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                      const WeightingFunctor<Dimension> weighting,
                      const int maxIterations,
                      const double tolerance) {

  typedef Mesh<Dimension> MeshType;
  typedef typename MeshType::Zone Zone;
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Grab the state.
  FieldList<Dimension, Vector> position = dataBase.globalPosition();
  FieldList<Dimension, SymTensor> H = dataBase.globalHfield();
  FieldList<Dimension, Scalar> weightSum = dataBase.newGlobalFieldList(0.0, "weightSum");
  FieldList<Dimension, Vector> delta = dataBase.newGlobalFieldList(Vector::zero, "delta");
  const Vector& xmin = boundary.xmin();
  const Vector& xmax = boundary.xmax();
  const double stopTol = tolerance*((xmax - xmin).maxAbsElement());
  const unsigned numNodeLists = dataBase.numNodeLists();

  // Iterate until we either hit the maximum number of iterations or hit the convergence
  // tolerance.
  double maxDelta = 2.0*stopTol;
  int iter = 0;
  while (iter < maxIterations and maxDelta > stopTol) {
    ++iter;
    weightSum = 0.0;
    delta = Vector::zero;
    maxDelta = 0.0;

    // Read out the node positions to a flat list.
    vector<Vector> generators;
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      copy(position[nodeListi]->begin(), 
           position[nodeListi]->begin() + n,
           back_inserter(generators));
    }

    // Generate the tessellation.
    MeshType mesh(generators, boundary);

    // Apply centroidal filtering.
    unsigned k = 0;
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        const Zone& zonei = mesh.zone(k);
        const vector<unsigned>& nodes = zonei.nodeIDs();
        ++k;

        // Compute a weighted centroid.
        double Wsum = 0.0;
        Vector centroid;
        for (vector<unsigned>::const_iterator itr = nodes.begin();
             itr != nodes.end();
             ++itr) {
          const Vector rj = mesh.node(*itr).position();
          const double Wj = weighting(rj, boundary);
          Wsum += Wj;
          centroid += Wj*rj;
        }
        CHECK(Wsum > 0.0);
        centroid /= Wsum;
        centroid = zonei.position();

        // Compute the delta to move toward the centroid.
        const Vector& ri = position(nodeListi, i);
        delta(nodeListi, i) = 0.5*(ri + centroid) - ri;

        // If we moved outside the boundary, go back!
        if (not boundary.contains(ri + delta(nodeListi, i))) {
          const Vector p = boundary.closestPoint(ri);
          delta(nodeListi, i) = 0.5*(ri + p) - ri;
        }
        maxDelta = std::max(maxDelta, delta(nodeListi, i).magnitude());
      }
    }
    if (Process::getRank() == 0) cerr << "relaxNodeDistribution iteration " << iter << " maxDelta=" << maxDelta << endl;

    // Apply the change.
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        position(nodeListi, i) += delta(nodeListi, i);
        CHECK(boundary.contains(position(nodeListi, i)));
      }
    }
  }
}

}

