//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/Timer.hh"
#include "Geometry/polyclipper.hh"

extern Timer TIME_computeVoronoiVolume2d;

#include <algorithm>
#include <utility>
#include <ctime>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// += for vector<vector<Plane>>
//
// Necessary to support reduction operations with
// FieldList<vector<vector<Plane>>> as a very specialized thing in this case.
//------------------------------------------------------------------------------
template<typename T>
inline
std::vector<std::vector<T>>&
operator+=(std::vector<std::vector<T>>& a, const std::vector<std::vector<T>>& b) {
  CHECK(a.size() == b.size());
  for (auto i = 0; i < a.size(); ++i) a[i].insert(a[i].end(), b[i].begin(), b[i].end());
  return a;
}

namespace Spheral {

using namespace FastMath;

namespace {  // anonymous namespace
// //------------------------------------------------------------------------------
// // A special comparator to sort r2d planes by distance.
// //------------------------------------------------------------------------------
// inline
// bool compareR2Dplanes(const r2d_plane& lhs, const r2d_plane& rhs) {
//   return lhs.d < rhs.d;
// }

// //------------------------------------------------------------------------------
// // Find the 1D extent of an R2D cell along the given direction.
// //------------------------------------------------------------------------------
// inline
// void findPolygonExtent(double& xmin, double& xmax, const Dim<2>::Vector& nhat, const r2d_poly& celli) {
//   REQUIRE(fuzzyEqual(nhat.magnitude(), 1.0));
//   double xi;
//   xmin = 0.0;
//   xmax = 0.0;
//   for (unsigned i = 0; i != celli.nverts; ++i) {
//     xi = (celli.verts[i].pos.x * nhat.x() +
//           celli.verts[i].pos.y * nhat.y());
//     xmin = std::min(xmin, xi);
//     xmax = std::max(xmax, xi);
//   }
//   xmin = std::min(0.0, xmin);
//   xmax = std::max(0.0, xmax);
// }

// // This version fits an ellipse and slices in the chosen direction.
// void findPolygonExtent(double& xmin, double& xmax, const Dim<2>::Vector& nhat, r2d_poly& celli) {
//   REQUIRE(fuzzyEqual(nhat.magnitude(), 1.0));
//   r2d_real moms[6];
//   r2d_reduce(&celli, moms, 2);
//   const Dim<2>::SymTensor G2(moms[3], moms[4], moms[4], moms[5]);
//   Dim<2>::SymTensor G = G2.sqrt();
//   // G *= sqrt(moms[0]/G.Determinant());
//   xmax = (G*nhat).magnitude();
//   xmin = -xmax;
// }

// //------------------------------------------------------------------------------
// // Worker function to clip by neighors in the given NodeList.
// // Returns whether the geometry was affected by this operation.
// //------------------------------------------------------------------------------
// bool clipByNeighbors(r2d_poly& celli,
//                      FieldList<Dim<2>, int>& surfacePoint,
//                      FieldList<Dim<2>, vector<Dim<2>::Vector>>& etaVoidPoints,
//                      const bool haveWeights,
//                      const bool returnSurface,
//                      const Dim<2>::Scalar rin,
//                      const vector<vector<int>>& fullConnectivity,
//                      const FieldList<Dim<2>, Dim<2>::Vector>& position,
//                      const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
//                      const unsigned nodeListi, 
//                      const unsigned i,
//                      const unsigned nodeListj) {

//   typedef Dim<2>::Scalar Scalar;
//   typedef Dim<2>::Vector Vector;
//   typedef Dim<2>::SymTensor SymTensor;
//   typedef Dim<2>::FacetedVolume FacetedVolume;
//   typedef Dim<2>::FacetedVolume::Facet Facet;

//   // Get the starting volume.
//   r2d_real vol0, vol1, vol2;
//   r2d_reduce(&celli, &vol0, 0);
//   CHECK(vol0 > 0.0);

//   // Check for multimaterial.
//   if (returnSurface and nodeListj != nodeListi and not fullConnectivity[nodeListj].empty()) surfacePoint(nodeListi, i) |= (1 << (nodeListj + 1));

//   // Build the clipping planes.
//   const auto& ri = position(nodeListi, i);
//   const auto  weighti = haveWeights ? weight(nodeListi, i) : 1.0;
//   vector<r2d_plane> pairPlanes, voidPlanes;
//   for (auto jItr = fullConnectivity[nodeListj].begin();
//        jItr != fullConnectivity[nodeListj].end();
//        ++jItr) {
//     const auto  j = *jItr;
//     const auto& rj = position(nodeListj, j);
//     const auto  weightj = haveWeights ? weight(nodeListj, j) : 1.0;

//     // Build the planes for the clipping half-spaces.
//     const auto rij = ri - rj;
//     const auto nhat = rij.unitVector();
//     const auto wij = weighti/(weightj + weighti);
//     if (voidPoint(nodeListj, j) == 0) {
//       pairPlanes.push_back(r2d_plane());
//       pairPlanes.back().n.x = nhat.x();
//       pairPlanes.back().n.y = nhat.y();
//       pairPlanes.back().d = wij*rij.magnitude();
//     } else {
//       voidPlanes.push_back(r2d_plane());
//       voidPlanes.back().n.x = nhat.x();
//       voidPlanes.back().n.y = nhat.y();
//       voidPlanes.back().d = wij*rij.magnitude();
//     }
//   }

//   // Sort the planes by distance -- lets us clip more efficiently.
//   std::sort(pairPlanes.begin(), pairPlanes.end(), compareR2Dplanes);
//   std::sort(voidPlanes.begin(), voidPlanes.end(), compareR2Dplanes);

//   // Clip the local cell.
//   r2d_clip(&celli, &pairPlanes[0], pairPlanes.size());
//   CHECK(celli.nverts > 0);
//   r2d_reduce(&celli, &vol1, 0);
//   CHECK(vol1 > 0.0);

//   // If there are void planes, check if they do anything.
//   if (not voidPlanes.empty()) {
//     const auto nvoid = voidPlanes.size();
//     for (auto k = 0; k != nvoid; ++k) {
//       r2d_clip(&celli, &voidPlanes[k], 1);
//       CHECK(celli.nverts > 0);
//       r2d_reduce(&celli, &vol2, 0);
//       if (vol2 < vol1) {
//         vol1 = vol2;
//         if (returnSurface) {
//           surfacePoint(nodeListi, i) |= 1;
//           etaVoidPoints(nodeListi, i).push_back(-0.5*rin*Vector(voidPlanes[k].n.x,
//                                                                 voidPlanes[k].n.y));
//         }
//       }
//     }
//   }

//   // Return whether we actually affected the cell.
//   return vol1 < vol0;
// }

//------------------------------------------------------------------------------
// Find the 1D extent of a polygon along the given direction.
//------------------------------------------------------------------------------
inline
void findPolygonExtent(double& xmin, double& xmax, 
                       const Dim<2>::Vector& nhat, 
                       const PolyClipper::Polygon& celli) {
  REQUIRE(fuzzyEqual(nhat.magnitude(), 1.0));
  double xi;
  xmin = 0.0;
  xmax = 0.0;
  for (const auto& vi: celli) {
    xi = vi.position.dot(nhat);
    xmin = std::min(xmin, xi);
    xmax = std::max(xmax, xi);
  }
  xmin = std::min(0.0, xmin);
  xmax = std::max(0.0, xmax);
}

}           // anonymous namespace

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
void
computeVoronoiVolume(const FieldList<Dim<2>, Dim<2>::Vector>& position,
                     const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                     const ConnectivityMap<Dim<2> >& connectivityMap,
                     const FieldList<Dim<2>, Dim<2>::SymTensor>& damage,
                     const std::vector<Dim<2>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<2>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<2>>*>& boundaries,
                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                     FieldList<Dim<2>, int>& surfacePoint,
                     FieldList<Dim<2>, Dim<2>::Scalar>& vol,
                     FieldList<Dim<2>, Dim<2>::Vector>& deltaMedian,
                     FieldList<Dim<2>, vector<Dim<2>::Vector>>& etaVoidPoints,
                     FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells,
                     FieldList<Dim<2>, std::vector<CellFaceFlag>>& cellFaceFlags) {

  TIME_computeVoronoiVolume2d.start();

  // Pre-conditions
  REQUIRE(facetedBoundaries.size() == 0 or facetedBoundaries.size() == position.size());
  REQUIRE(vol.size() == position.size());
  REQUIRE(deltaMedian.size() == position.size());
  REQUIRE(holes.size() == facetedBoundaries.size());

  typedef Dim<2> Dimension;
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef PolyClipper::Plane2d Plane;

  const auto numGens = position.numNodes();
  const auto numNodeLists = position.size();
  const auto numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());
  const auto haveFacetedBoundaries = facetedBoundaries.size() == numNodeLists;
  const auto haveBoundaries = not boundaries.empty();
  const auto haveWeights = weight.size() == numNodeLists;
  const auto haveDamage = damage.size() == numNodeLists;
  const auto returnSurface = surfacePoint.size() == numNodeLists;
  const auto returnCells = cells.size() == numNodeLists;

  if (returnSurface) {
    surfacePoint = 0;
    etaVoidPoints = vector<Vector>();
  }

  if (numGensGlobal > 0) {

    // The criterion for the initial starting sphere of a point in eta space.
    const auto rin = 2.0/vol[0]->nodeListPtr()->nodesPerSmoothingScale();

    // Unit circle as template shape.
    const auto nverts = 18;
    PolyClipper::Polygon cell0;
    {
      const auto dtheta = 2.0*M_PI/nverts;
      vector<Vector> verts0(nverts);
      vector<vector<unsigned>> facets0(nverts, vector<unsigned>(2));
      for (auto j = 0; j != nverts; ++j) {
        const auto theta = j*dtheta;
        verts0[j].x(cos(theta));
        verts0[j].y(sin(theta));
        facets0[j][0] = j;
        facets0[j][1] = (j + 1) % nverts;
      }
      PolyClipper::convertToPolygon(cell0, FacetedVolume(verts0, facets0));
    }

    // We'll need to hang onto the PolyClipper cells and any per cell void points.
    FieldList<Dim<2>, PolyClipper::Polygon> polycells(FieldStorageType::CopyFields);
    FieldList<Dim<2>, vector<vector<Plane>>> pairPlanes(FieldStorageType::CopyFields);
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      polycells.appendNewField("polycells", vol[nodeListi]->nodeList(), cell0);
      pairPlanes.appendNewField("pair planes", vol[nodeListi]->nodeList(), vector<vector<Plane>>(numNodeLists));
    }

    //==========================================================================
    // First pass: clip by any faceted boundaries/holes.
    for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
      const auto ni = polycells[nodeListi]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0; i < ni; ++i) {

        // Initialize the per point polygon by scaling by its H extent
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hinv = Hi.Inverse();
        for (auto& v: polycells(nodeListi, i)) v.position = 1.1*rin*Hinv*v.position;

        // Clip by any faceted boundaries first.
        if (haveFacetedBoundaries) {
          vector<Plane> boundPlanes;
          const auto& facets = facetedBoundaries[nodeListi].facets();
          CHECK(facetedBoundaries[nodeListi].contains(ri, false));
          for (const auto& facet: facets) {
            const auto p = facet.closestPoint(ri);
            auto rji = p - ri;
            if ((Hi*rji).magnitude2() < rin*rin) {
              Vector nhat;
              if (rji.magnitude() < 1.0e-5*facet.area()) {
                rji.Zero();
                nhat = -facet.normal();
              } else {
                nhat = -rji.unitVector();
              }
              boundPlanes.push_back(Plane(rji, nhat));
            }
          }

          // Same thing with holes.
          for (const auto& hole: holes[nodeListi]) {
            CHECK(not hole.contains(ri, false));
            const auto& facets = hole.facets();
            for (const auto& facet: facets) {
              const auto p = facet.closestPoint(ri);
              auto rji = p - ri;
              if ((Hi*rji).magnitude2() < rin*rin) {
                Vector nhat;
                if (rji.magnitude2() < 1.0e-5*facet.area()) {
                  rji.Zero();
                  nhat = facet.normal();
                } else {
                  nhat = -rji.unitVector();
                }
                boundPlanes.push_back(Plane(rji, nhat));
              }
            }
          }

          // Sort the planes by distance -- lets us clip more efficiently.
          std::sort(boundPlanes.begin(), boundPlanes.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });

          // Clip by the planes thus far.
          PolyClipper::clipPolygon(polycells(nodeListi, i), boundPlanes);
        }
      }
    }

    //==========================================================================
    // Second pass: clip by neighbor points.  Note we have to keep track of
    // which NodeLists actually clip each polygon in order to detect material
    // surfaces.
    const auto& pairs = connectivityMap.nodePairList();
    const auto  npairs = pairs.size();
#pragma omp parallel
    {
      // Thread private scratch variables
      int i, j, nodeListi, nodeListj;
      auto pairPlanes_thread = pairPlanes.threadCopy();

#pragma omp for
      for (auto kk = 0; kk < npairs; ++kk) {
        i = pairs[kk].i_node;
        j = pairs[kk].j_node;
        nodeListi = pairs[kk].i_list;
        nodeListj = pairs[kk].j_list;

        // State of node i
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hinv = Hi.Inverse();
        const auto  weighti = haveWeights ? weight(nodeListi, i) : 1.0;
        auto&       pairPlanesi = pairPlanes_thread(nodeListi, i);

        // State of node j
        const auto& rj = position(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  Hjnv = Hj.Inverse();
        const auto  weightj = haveWeights ? weight(nodeListj, j) : 1.0;
        auto&       pairPlanesj = pairPlanes_thread(nodeListj, j);

        // Build the planes for the clipping half-spaces.
        const auto rji = rj - ri;
        const auto nhat = -rji.unitVector();
        const auto wij = weighti/(weighti + weightj);
        const auto wji = weightj/(weighti + weightj);
        pairPlanesi[nodeListj].push_back(Plane( wij*rji,  nhat));
        pairPlanesj[nodeListi].push_back(Plane(-wji*rji, -nhat));
      }

      // Collect the pair planes across threads.
      pairPlanes_thread.threadReduce();
    } // OMP parallel

    // Clip by the neighbors.
    for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
      const auto ni = polycells[nodeListi]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0; i < ni; ++i) {
        const auto& Hi = H(nodeListi, i);
        const auto  Hinvi = Hi.Inverse();
        auto& celli = polycells(nodeListi, i);
        auto& pairPlanesi = pairPlanes(nodeListi, i);
        if (returnSurface) {

          // If we're returning the surface we have to pay attention to which materials actually clip this point.
          double vol0, vol1;
          Vector cent;
          PolyClipper::moments(vol0, cent, celli);
          for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
            std::sort(pairPlanesi[nodeListj].begin(), pairPlanesi[nodeListj].end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
            PolyClipper::clipPolygon(celli, pairPlanesi[nodeListj]);
            PolyClipper::moments(vol1, cent, celli);
            if (vol1 < vol0) {
              vol0 = vol1;
              if (nodeListj != nodeListi) surfacePoint(nodeListi, i) |= (1 << (nodeListj + 1));
            }
          }

        } else {

          // Otherwise just clip by all the pair planes at once.
          vector<Plane> planes = pairPlanesi[0];
          for (auto nodeListj = 1; nodeListj < numNodeLists; ++nodeListj) planes += pairPlanesi[nodeListj];
          std::sort(planes.begin(), planes.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
          PolyClipper::clipPolygon(celli, planes);

        }

        // Check if the final polygon is entirely within our "interior" check radius.  Otherwise,
        // time to make void points.
        auto nvoid = 0;
        Vector etaVoidAvg;
        for (const auto& vert: celli) {
          const auto peta = Hi*vert.position;
          if (peta.magnitude2() > rin*rin) {
            ++nvoid;
            etaVoidAvg += peta.unitVector();
          }
        }

        // Is this point on a free surface?
        // If so, we have to generate void points.
        if (nvoid > 0) {

          // Reduce the number of void points for this point to a reasonable stencil -- up to 4 for 2D.
          CHECK(etaVoidPoints(nodeListi, i).empty());
          const auto thetaVoidAvg = atan2(etaVoidAvg.y(), etaVoidAvg.x());
          const auto nv = max(1U, min(4U, unsigned(4.0*float(nvoid)/float(nverts))));
          for (unsigned k = 0; k != nv; ++k) {
            const auto theta = thetaVoidAvg + (0.5*k - 0.25*(nv - 1))*M_PI;
            etaVoidPoints(nodeListi, i).push_back(Vector(0.5*rin*cos(theta), 0.5*rin*sin(theta)));
          }
          CHECK(etaVoidPoints(nodeListi, i).size() == nv);
        }

        // If this point is sufficiently damaged, we also create void points along the damaged directions.
        if (haveDamage and damage(nodeListi, i).Trace() > 1.0 - 1.0e-5) {
          const auto ev = damage(nodeListi, i).eigenVectors();
          for (auto jdim = 0; jdim < Dimension::nDim; ++jdim) {
            if (ev.eigenValues(jdim) > 1.0 - 1.0e-5) {
              const auto evecj = ev.eigenVectors.getColumn(jdim);
              etaVoidPoints(nodeListi, i).push_back(-0.5*rin*evecj);
              etaVoidPoints(nodeListi, i).push_back( 0.5*rin*evecj);
            }
          }
          surfacePoint(nodeListi, i) |= 1;

          // We reuse the first element of pairPlanes to build up the void planes.
          // Initialize with each points own void point planes.  We'll add the void planes from neighbors
          // in the next phase.
          pairPlanes(nodeListi, i) = vector<vector<Plane>>(1);
          for (const auto& etaVoid: etaVoidPoints(nodeListi, i)) {
            const auto rji = Hinvi*etaVoid;
            const auto nhat = -rji.unitVector();
            pairPlanes(nodeListi, i)[0].push_back(Plane(rji, nhat));
          }
        }
      }  // end over i (OMP as well)
    }    // end over NodeLists

    // Apply boundary conditions to the void points.
    if (not boundaries.empty()) {
      for (const auto bc: boundaries) bc->applyFieldListGhostBoundary(etaVoidPoints);
      for (const auto bc: boundaries) bc->finalizeGhostBoundary();
    }

    //==========================================================================
    // Third pass: clip by any void points.
    // Add the void points from neighbors.
#pragma omp parallel
    {
      // Thread private scratch variables
      int i, j, nodeListi, nodeListj;
      auto pairPlanes_thread = pairPlanes.threadCopy();

#pragma omp for
      for (auto kk = 0; kk < npairs; ++kk) {
        i = pairs[kk].i_node;
        j = pairs[kk].j_node;
        nodeListi = pairs[kk].i_list;
        nodeListj = pairs[kk].j_list;

        // State of node i
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hinvi = Hi.Inverse();
        auto&       pairPlanesi = pairPlanes_thread(nodeListi, i)[0];

        // State of node j
        const auto& rj = position(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  Hinvj = Hj.Inverse();
        auto&       pairPlanesj = pairPlanes_thread(nodeListj, j)[0];

        for (const auto& etaVoid: etaVoidPoints(nodeListj, j)) {
          const auto rji = Hinvj*etaVoid + 0.5*(rj - ri);
          const auto nhat = -rji.unitVector();
          pairPlanesi.push_back(Plane(rji, nhat));
        }

        for (const auto& etaVoid: etaVoidPoints(nodeListi, i)) {
          const auto rij = Hinvi*etaVoid + 0.5*(ri - rj);
          const auto nhat = -rij.unitVector();
          pairPlanesj.push_back(Plane(rij, nhat));
        }
      }

      pairPlanes_thread.threadReduce();

    } // OMP parallel

    // Now we can do the void clipping, compute the final volumes, and finalize
    // surface detection.
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = vol[nodeListi]->numInternalElements();
#pragma omp parallel for
      for (auto i = 0; i < n; ++i) {
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hinvi = Hi.Inverse();
        auto&       voidPlanesi = pairPlanes(nodeListi, i)[0];
        auto&       celli = polycells(nodeListi, i);

        // Clip by the void planes, compute the volume, and check if the void surface had any effect.
        double vol0, vol1;
        PolyClipper::moments(vol0, deltaMedian(nodeListi, i), celli);
        std::sort(voidPlanesi.begin(), voidPlanesi.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
        PolyClipper::clipPolygon(celli, voidPlanesi);
        CHECK(not celli.empty());
        PolyClipper::moments(vol1, deltaMedian(nodeListi, i), celli);

        // We only use the volume result if interior.
        const bool interior = (vol1 >= vol0);
        if (interior) {
          vol(nodeListi, i) = vol1;
        } else if (returnSurface) {
          surfacePoint(nodeListi, i) |= 1;
        }

        // Check if the candidate motion is still in the boundary.  If not, project back.
        if (haveFacetedBoundaries) {
          if (not facetedBoundaries[nodeListi].contains(ri + deltaMedian(nodeListi, i), false)) {
            deltaMedian(nodeListi, i) = facetedBoundaries[nodeListi].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
          }
          for (unsigned ihole = 0; ihole != holes[nodeListi].size(); ++ihole) {
            if (holes[nodeListi][ihole].contains(ri + deltaMedian(nodeListi, i), false)) {
              deltaMedian(nodeListi, i) = holes[nodeListi][ihole].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
            }
          }
        }

        // If requested, we can return the cell geometries.
        if (returnCells) {
          PolyClipper::collapseDegenerates(celli, 1.0e-10);
          PolyClipper::convertFromPolygon(cells(nodeListi, i), celli);
          cells(nodeListi, i) += ri;
        }
      }
    }
  }

  TIME_computeVoronoiVolume2d.stop();
}
    
}
