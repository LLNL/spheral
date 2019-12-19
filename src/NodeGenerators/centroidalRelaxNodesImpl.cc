//------------------------------------------------------------------------------
// Implement Lloyd's algorithm for centroidal relaxation of fluid points.
//------------------------------------------------------------------------------
#include "centroidalRelaxNodesImpl.hh"
#include "CRKSPH/computeVoronoiVolume.hh"
#include "CRKSPH/computeCRKSPHMoments.hh"
#include "CRKSPH/computeCRKSPHCorrections.hh"
#include "CRKSPH/gradientCRKSPH.hh"
#include "NodeList/ASPHSmoothingScale.hh"
#include "Utilities/iterateIdealH.hh"

#include <ctime>
using std::vector;
using std::tuple;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

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
                         const double fracTol,
                         const CRKOrder correctionOrder,
                         const double centroidFrac,
                         FieldList<Dimension, double>& vol,
                         FieldList<Dimension, int>& surfacePoint,
                         FieldList<Dimension, typename Dimension::FacetedVolume>& cells) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Read some state.
  const auto numNodeLists = db.numNodeLists();
  const bool useBounds = (volumeBoundaries.size() == numNodeLists);
  const bool useHoles = (holes.size() == numNodeLists);
  auto pos = db.fluidPosition();
  auto H = db.fluidHfield();
  auto mass = db.fluidMass();
  auto rhof = db.fluidMassDensity();
  auto D = db.solidEffectiveDamage();

  // Prepare the storage for the point-wise fields.
  auto gradRhof = db.newFluidFieldList(Vector::zero, "mass density gradient");
  auto deltaCentroid = db.newFluidFieldList(Vector::zero, "delta centroid");
  auto A = db.newFluidFieldList(0.0, "A");
  auto B = db.newFluidFieldList(Vector::zero, "B");
  auto C = db.newFluidFieldList(Tensor::zero, "B");
  auto gradA = db.newFluidFieldList(Vector::zero, "gradA");
  auto gradB = db.newFluidFieldList(Tensor::zero, "gradB");
  auto gradC = db.newFluidFieldList(ThirdRankTensor::zero, "gradC");
  auto m0 = db.newFluidFieldList(0.0, "m0");
  auto m1 = db.newFluidFieldList(Vector::zero, "m1");
  auto m2 = db.newFluidFieldList(SymTensor::zero, "m2");
  auto m3 = db.newFluidFieldList(ThirdRankTensor::zero, "m3");
  auto m4 = db.newFluidFieldList(FourthRankTensor::zero, "m4");
  auto gradm0 = db.newFluidFieldList(Vector::zero, "gradm0");
  auto gradm1 = db.newFluidFieldList(Tensor::zero, "gradm1");
  auto gradm2 = db.newFluidFieldList(ThirdRankTensor::zero, "gradm2");
  auto gradm3 = db.newFluidFieldList(FourthRankTensor::zero, "gradm3");
  auto gradm4 = db.newFluidFieldList(FifthRankTensor::zero, "gradm4");

  // Temporary until we decide to propagate void info to this method.
  auto etaVoidPoints = db.newFluidFieldList(vector<Vector>(), "eta void points");

  // Make a dummy set of cells so we don't ask computeVoronoiVolume to compute the return FacetedVolumes every step.
  FieldList<Dimension, FacetedVolume> dummyCells;
  FieldList<Dimension, vector<CellFaceFlag>> cellFaceFlags;

  // Make sure the density starts out consistently, and kick-start the volume using m/rho.
  for (auto nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
    const auto n = rhof[nodeListi]->numInternalElements();
    for (auto i = 0U; i != n; ++i) {
      rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
      CHECK(rhof(nodeListi, i) > 0.0);
      vol(nodeListi, i) = mass(nodeListi, i)/rhof(nodeListi, i);
    }
  }

  // // Update the H tensors a bit.
  // iterateIdealH(db, boundaries, W, ASPHSmoothingScale<Dimension>(), 5);

  // Iterate until we converge or max out.
  unsigned iter = 0;
  double avgdelta = 2.0*fracTol;
  while (iter < 2 or (iter < maxIterations and (avgdelta > fracTol))) { //  or mass.min() == 0.0))) {
    std::clock_t t0 = std::clock();
    iter += 1;

    // Remove any old ghost nodes info, and update the mass density
    {
      unsigned nodeListi = 0U;
      for (auto& nodes: db.fluidNodeListPtrs()) {
        nodes->numGhostNodes(0);
        nodes->neighbor().updateNodes();
        if (not rhoConst) for (unsigned i = 0; i != nodes->numInternalNodes(); ++i) rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
      }
      ++nodeListi;
    }

    // Create the new ghost nodes.
    for (auto& bc: boundaries) bc->setAllGhostNodes(db);
    for (auto& bc: boundaries) bc->finalizeGhostBoundary();
    for (auto& nodes: db.fluidNodeListPtrs()) nodes->neighbor().updateNodes();

    // Compute the new connectivity.
    std::clock_t tcm = std::clock();
    db.updateConnectivityMap(false, false);
    const auto& cm = db.connectivityMap();
    tcm = std::clock() - tcm;
    { // BLAGO
      double avgneighbors = 0.0;
      unsigned ntot = 0;
      for (auto nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
        const auto n = rhof[nodeListi]->numInternalElements();
        ntot += n;
        for (auto i = 0U; i != n; ++i) avgneighbors += cm.numNeighborsForNode(nodeListi, i);
      }
      ntot = allReduce(ntot, MPI_SUM, Communicator::communicator());
      avgneighbors = allReduce(avgneighbors, MPI_SUM, Communicator::communicator())/ntot;
      if (Process::getRank() == 0) cerr << "Avergage number of neighbors per node: " << avgneighbors << " " << ntot << endl;
    } // BLAGO

    // Compute the new volumes and centroids (note this uses the old rho gradient, not quite right,
    // but expedient/efficient).
    std::clock_t tvoro = std::clock();
    computeVoronoiVolume(pos, H, cm, D, volumeBoundaries, holes, boundaries,
                         FieldList<Dimension, typename Dimension::Scalar>(),  // no weights
                         surfacePoint, vol, deltaCentroid, etaVoidPoints, dummyCells, cellFaceFlags);
    tvoro = std::clock() - tvoro;
     
    // Apply boundary conditions.
    for (auto& bc: boundaries) {
      bc->applyFieldListGhostBoundary(vol);
      bc->applyFieldListGhostBoundary(rhof);
    }
    for (auto& bc: boundaries) bc->finalizeGhostBoundary();

    // If the density is constant we can entirely skip finding the gradient.
    if (not rhoConst) {

      // If the user provided a gradrho method, we can use it.  Otherwise we need to numerically evaluate
      // the density gradient.
      if (useGradRhoFunc) {
        for (unsigned nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
          const auto n = rhof[nodeListi]->numInternalElements();
          for (auto i = 0U; i != n; ++i) {
            gradRhof(nodeListi, i) = gradrhofunc(pos(nodeListi, i));
          }
        }

      } else {
        // Use RK to numerically compute the new mass density gradient.
        computeCRKSPHMoments(cm, W, vol, pos, H, correctionOrder, NodeCoupling(),
                                          m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
        computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, H, surfacePoint, correctionOrder,
                                              A, B, C, gradA, gradB, gradC);
        gradRhof.assignFields(gradientCRKSPH(rhof, pos, vol, H, A, B, C, gradA, gradB, gradC, cm, correctionOrder, W));
      }
    }
     
    // Displace the points and update point masses.
    avgdelta = 0.0;
    for (unsigned nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = rhof[nodeListi]->numInternalElements();
      const auto hmax = rhof[nodeListi]->nodeListPtr()->hmax();
      for (auto i = 0U; i != n; ++i) {
        auto delta = centroidFrac * deltaCentroid(nodeListi, i);
        if (useBounds) {
          while (not volumeBoundaries[nodeListi].contains(pos(nodeListi, i) + delta, false)) delta *= 0.9;
          if (useHoles) {
            for (const auto& hole: holes[nodeListi]) {
              while (hole.contains(pos(nodeListi, i) + delta, false)) delta *= 0.9;
            }
          }
        }
        if (vol(nodeListi, i) > 0.0) {
          avgdelta += delta.magnitude()/Dimension::rootnu(vol(nodeListi, i));
          H(nodeListi, i) = SymTensor::one / std::min(hmax, 2.0*Dimension::rootnu(vol(nodeListi, i)));  // Not correct, but hopefully good enough for our iterative Voronoi purposes.
        }
        pos(nodeListi, i) += delta;
        rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
        mass(nodeListi, i) = rhof(nodeListi,i)*vol(nodeListi,i);
      }
    }
    avgdelta = (allReduce(avgdelta, MPI_SUM, Communicator::communicator())/
                allReduce(db.numInternalNodes(), MPI_SUM, Communicator::communicator()));
    if (Process::getRank() == 0) cerr << "centroidalRelaxNodes iteration " << iter << ", avg delta frac " << avgdelta 
                                      << ", required " << ((std::clock() - t0)/CLOCKS_PER_SEC) 
                                      << " seconds (" << (tvoro/CLOCKS_PER_SEC) << " in computeVoronoiVolume, "
                                      << (tcm/CLOCKS_PER_SEC) << " in ConnectivityMap)." << endl;
        
    // // Update the H tensors a bit.
    // iterateIdealH(db, boundaries, W, ASPHSmoothingScale<Dimension>(), 2);
  }

  // If requested to return the FacetedVolumes, make one last call to fill 'em in.
  if (cells.size() > 0) {
    const auto& cm = db.connectivityMap();
    computeVoronoiVolume(pos, H, cm, D, volumeBoundaries, holes, boundaries,
                         FieldList<Dimension, typename Dimension::Scalar>(),  // no weights
                         surfacePoint, vol, deltaCentroid, etaVoidPoints, cells, cellFaceFlags);
  }

  // Return how many iterations we actually took.
  return iter;
}

}
