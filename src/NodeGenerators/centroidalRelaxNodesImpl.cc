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

namespace Spheral {

using namespace std;
using std::min;
using std::max;
using std::abs;

using KernelSpace::TableKernel;
using BoundarySpace::Boundary;
using DataBaseSpace::DataBase;

template<typename Dimension>
unsigned
centroidalRelaxNodesImpl(DataBaseSpace::DataBase<Dimension>& db,
                         const std::vector<typename Dimension::FacetedVolume>& volumeBoundaries,
                         const std::vector<std::vector<typename Dimension::FacetedVolume> >& holes,
                         const KernelSpace::TableKernel<Dimension>& W,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, double>& rhofunc,
                         const PythonBoundFunctors::SpheralFunctor<typename Dimension::Vector, typename Dimension::Vector>& gradrhofunc,
                         const bool rhoConst,
                         const bool useGradRhoFunc,
                         std::vector<BoundarySpace::Boundary<Dimension>*>& boundaries,
                         const unsigned maxIterations,
                         const double fracTol,
                         const CRKSPHSpace::CRKOrder correctionOrder,
                         const double centroidFrac,
                         FieldSpace::FieldList<Dimension, double>& vol,
                         FieldSpace::FieldList<Dimension, int>& surfacePoint,
                         FieldSpace::FieldList<Dimension, typename Dimension::FacetedVolume>& cells) {

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

  // Make a dummy set of cells so we don't ask computeVoronoiVolume to compute the return FacetedVolumes every step.
  FieldList<Dimension, FacetedVolume> dummyCells;

  // Make sure the density starts out consistently, and kick-start the volume using m/rho.
  for (auto nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
    const auto n = rhof[nodeListi]->numInternalElements();
    for (auto i = 0U; i != n; ++i) {
      rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
      CHECK(rhof(nodeListi, i) > 0.0);
      vol(nodeListi, i) = mass(nodeListi, i)/rhof(nodeListi, i);
    }
  }

  // Iterate until we converge or max out.
  unsigned iter = 0;
  double avgdelta = 2.0*fracTol;
  while (iter < 2 or (iter < maxIterations and (avgdelta > fracTol))) { //  or mass.min() == 0.0))) {
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
    db.updateConnectivityMap(false);
    const auto& cm = db.connectivityMap();

    // Compute the new volumes and centroids (note this uses the old rho gradient, not quite right,
    // but expedient/efficient).
    CRKSPHSpace::computeVoronoiVolume(pos, H, rhof, gradRhof, cm, W.kernelExtent(), volumeBoundaries, holes, 
                                      FieldList<Dimension, typename Dimension::Scalar>(),  // no weights
                                      surfacePoint, vol, deltaCentroid, dummyCells);
     
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
        CRKSPHSpace::computeCRKSPHMoments(cm, W, vol, pos, H, correctionOrder, NodeCoupling(),
                                          m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
        CRKSPHSpace::computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, H, correctionOrder,
                                              A, B, C, gradA, gradB, gradC);
        gradRhof.assignFields(CRKSPHSpace::gradientCRKSPH(rhof, pos, vol, H, A, B, C, gradA, gradB, gradC, cm, correctionOrder, W));
      }
    }
     
    // Displace the points and update point masses.
    avgdelta = 0.0;
    for (unsigned nodeListi = 0U; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = rhof[nodeListi]->numInternalElements();
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
          H(nodeListi, i) = SymTensor::one / (2.0*Dimension::rootnu(vol(nodeListi, i)));  // Not correct, but hopefully good enough for our iterative Voronoi purposes.
        }
        pos(nodeListi, i) += delta;
        rhof(nodeListi, i) = rhofunc(pos(nodeListi, i));
        mass(nodeListi, i) = rhof(nodeListi,i)*vol(nodeListi,i);
      }
    }
    avgdelta = (allReduce(avgdelta, MPI_SUM, Communicator::communicator())/
                allReduce(db.numInternalNodes(), MPI_SUM, Communicator::communicator()));
    if (Process::getRank() == 0) cout << "centroidalRelaxNodes iteration " << iter << " avg delta frac " << avgdelta << endl;
        
    // // Update the H tensors a bit.
    // iterateIdealH(db, boundaries, W, NodeSpace::ASPHSmoothingScale<Dimension>(), 2);
  }

  // If requested to return the FacetedVolumes, make one last call to fill 'em in.
  if (cells.size() > 0) {
    const auto& cm = db.connectivityMap();
    CRKSPHSpace::computeVoronoiVolume(pos, H, rhof, gradRhof, cm, W.kernelExtent(), volumeBoundaries, holes, 
                                      FieldList<Dimension, typename Dimension::Scalar>(),  // no weights
                                      surfacePoint, vol, deltaCentroid, cells);
  }

  // Return how many iterations we actually took.
  return iter;
}

}
