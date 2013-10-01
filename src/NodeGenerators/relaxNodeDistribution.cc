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
                      const WeightingFunctor<Dimension> weighting,
                      const int maxIterations,
                      const double tolerance) {

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
  Vector rImage;
  Scalar wj, wjimage;
  while (iter < maxIterations and maxDelta > stopTol) {
    ++iter;
    weightSum = 0.0;
    delta = Vector::zero;
    maxDelta = 0.0;

    // Update the connectivity.
    for (typename DataBase<Dimension>::ConstNodeListIterator itr = dataBase.nodeListBegin();
         itr != dataBase.nodeListEnd();
         ++itr) (*itr)->neighbor().updateNodes();
    dataBase.updateConnectivityMap();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();

    // Walk the NodeLists.
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {

      // Iterate over the internal nodes in this NodeList.
      for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
           iItr != connectivityMap.end(nodeListi);
           ++iItr) {

        // Get the state for node i.
        const int i = *iItr;
        const Vector& ri = position(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

        // // Self-contribution.
        // weightSum(nodeListi, i) = wi*W.kernelValue(0.0, Hdeti);
        // delta(nodeListi, i) = Vector::zero;

        // // Check for any contribution from the boundary.
        // const Vector rj = boundary.closestPoint(ri);
        // const Scalar wj = weighting(rj, boundary);
        // const Vector rji = rj - ri;
        // const Scalar etai = (Hi*rji).magnitude();
        // const Scalar Wi = std::abs(W.gradValue(etai, 1.0));
        // delta(nodeListi, i) += wj*Wi * rji;
        // weightSum(nodeListi, i) += wj*Wi;

        // Walk the neighbors.
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

          // Connectivity of this node with this NodeList.  We only need to proceed if
          // there are some nodes in this list.
          const vector<int>& connectivity = fullConnectivity[nodeListj];
          if (connectivity.size() > 0) {
            for (vector<int>::const_iterator jItr = connectivity.begin();
                 jItr != connectivity.end();
                 ++jItr) {

              // State for node j.
              const int j = *jItr;
              const Vector& rj = position(nodeListj, j);
              weighting(rj, boundary, rImage, wj, wjimage);

              // Sum the contribution from this point.
              {
                const Vector rji = rj - ri;
                const Scalar etaMagi = (Hi*rji).magnitude();
                const Scalar etai = (Hi*rji).magnitude();
                const Scalar Wi = std::abs(W.gradValue(etai, 1.0));
                delta(nodeListi, i) += wj*Wi * rji;
                weightSum(nodeListi, i) += wj*Wi;
              }

              // Sum the contribution from the image.
              {
                const Vector rji = rImage - ri;
                const Scalar etaMagi = (Hi*rji).magnitude();
                const Scalar etai = (Hi*rji).magnitude();
                const Scalar Wi = std::abs(W.gradValue(etai, 1.0));
                delta(nodeListi, i) += wjimage*Wi * rji;
                weightSum(nodeListi, i) += wjimage*Wi;
              }
            }
          }
        }

        // Set the delta for this node.
        CHECK(weightSum(nodeListi, i) > 0.0);
        delta(nodeListi, i) /= weightSum(nodeListi, i);

        // If we moved outside the boundary, go back!
        if (not boundary.contains(ri + delta(nodeListi, i))) {
          const Vector p = boundary.closestPoint(ri);
          delta(nodeListi, i) = 0.5*(ri + p) - ri;
        }

        // Update the maximum delta.
        maxDelta = std::max(maxDelta, delta(nodeListi, i).magnitude());
      }
    }
    if (Process::getRank() == 0) cerr << "relaxNodeDistribution iteration " << iter << " maxDelta=" << maxDelta << endl;

    // Apply the change.
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
           iItr != connectivityMap.end(nodeListi);
           ++iItr) {
        const int i = *iItr;
        position(nodeListi, i) += delta(nodeListi, i);
        CHECK(boundary.contains(position(nodeListi, i)));
      }
    }

  }

}

}

