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
//                      const FieldList<Dim<2>, int>& voidPoint,
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
                     const FieldList<Dim<2>, Dim<2>::Scalar>& rho,
                     const FieldList<Dim<2>, Dim<2>::Vector>& gradRho,
                     const ConnectivityMap<Dim<2> >& connectivityMap,
                     const FieldList<Dim<2>, Dim<2>::SymTensor>& damage,
                     const std::vector<Dim<2>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<2>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<2>>*>& boundaries,
                     const FieldList<Dim<2>, Dim<2>::Scalar>& weight,
                     const FieldList<Dim<2>, int>& voidPoint,
                     FieldList<Dim<2>, int>& surfacePoint,
                     FieldList<Dim<2>, Dim<2>::Scalar>& vol,
                     FieldList<Dim<2>, Dim<2>::Vector>& deltaMedian,
                     FieldList<Dim<2>, vector<Dim<2>::Vector>>& etaVoidPoints,
                     FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells) {

  TIME_computeVoronoiVolume2d.start();

  // Pre-conditions
  REQUIRE(facetedBoundaries.size() == 0 or facetedBoundaries.size() == position.size());
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

    // We'll need to hang onto the PolyClipper cells.
    FieldList<Dim<2>, PolyClipper::Polygon> polycells(FieldStorageType::CopyFields);
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      polycells.appendNewField("polycells", vol[nodeListi]->nodeList(), PolyClipper::Polygon());
    }

    // First pass: clip by neighbors and generate void points.
    double vol0, voli;
    vector<Plane> pairPlanes;
    unsigned nvoid;
    Vector etaVoidAvg;
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = vol[nodeListi]->numInternalElements();
      const auto rin = 2.0/vol[nodeListi]->nodeListPtr()->nodesPerSmoothingScale();

#pragma omp parallel for                                        \
  private(pairPlanes, nvoid, etaVoidAvg, vol0, voli)
      for (auto i = 0; i < n; ++i) {
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hinv = Hi.Inverse();
        const auto  weighti = haveWeights ? weight(nodeListi, i) : 1.0;
        const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        auto&       celli = polycells(nodeListi, i);

        // Prepare to accumulate any void point positions.
        etaVoidAvg = Vector::zero;
        nvoid = 0;
        pairPlanes.clear();

        // Initialize our seed cell shape.
        celli = cell0;
        for (auto& v: celli) v.position = 1.1*rin*Hinv*v.position;

        // Clip by any faceted boundaries first.
        if (haveFacetedBoundaries) {
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
              pairPlanes.push_back(Plane(rji, nhat));
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
                pairPlanes.push_back(Plane(rji, nhat));
              }
            }
          }
        }

        // Add clipping planes from neighbors in our own NodeList.
        for (auto j: fullConnectivity[nodeListi]) {
          const auto& rj = position(nodeListi, j);
          const auto  weightj = haveWeights ? weight(nodeListi, j) : 1.0;

          // Build the planes for the clipping half-spaces.
          const auto rji = rj - ri;
          const auto nhat = -rji.unitVector();
          const auto wij = weighti/(weightj + weighti);
          pairPlanes.push_back(Plane(wij*rji, nhat));
        }

        // Sort the planes by distance -- lets us clip more efficiently.
        std::sort(pairPlanes.begin(), pairPlanes.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });

        // Clip by the planes thus far.
        PolyClipper::clipPolygon(celli, pairPlanes);
        CHECK(not celli.empty());
        PolyClipper::moments(voli, deltaMedian(nodeListi, i), celli);

        // Now clip by neighbors in other NodeLists.
        for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          if (nodeListj != nodeListi) {
            pairPlanes.clear();
            for (auto j: fullConnectivity[nodeListj]) {
              const auto& rj = position(nodeListj, j);
              const auto  weightj = haveWeights ? weight(nodeListj, j) : 1.0;
              const auto rji = rj - ri;
              const auto nhat = -rji.unitVector();
              const auto wij = weighti/(weightj + weighti);
              pairPlanes.push_back(Plane(wij*rji, nhat));
            }

            // Clip by this NodeList's planes, and see if the volume changed.
            std::sort(pairPlanes.begin(), pairPlanes.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
            PolyClipper::clipPolygon(celli, pairPlanes);
            if (returnSurface) {
              vol0 = voli;
              PolyClipper::moments(voli, deltaMedian(nodeListi, i), celli);
              if (voli < vol0) surfacePoint(nodeListi, i) |= (1 << (nodeListj + 1));
            }
          }
        }

        // Check if the final polygon is entirely within our "interior" check radius.
        bool interior = true;
        for (const auto& vert: celli) {
          const auto peta = Hi*vert.position;
          if (peta.magnitude2() > rin*rin) {
            interior = false;
            if (returnSurface) {
              surfacePoint(nodeListi, i) |= 1;
              ++nvoid;
              etaVoidAvg += peta.unitVector();
            }
          }
        }

        // Reduce the number of void points for this point to a reasonable stencil -- up to 4 for 2D.
        CHECK(etaVoidPoints(nodeListi, i).empty());
        if (returnSurface and not interior) {
          CHECK(nvoid > 0);
          const Scalar thetaVoidAvg = atan2(etaVoidAvg.y(), etaVoidAvg.x());
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
        }
      }
    }

    // Apply boundary conditions to the void points.
    if (not boundaries.empty()) {
      for (const auto bc: boundaries) bc->applyFieldListGhostBoundary(etaVoidPoints);
      for (const auto bc: boundaries) bc->finalizeGhostBoundary();
    }

    // Second pass: clip by any void points and compute the volumes.
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = vol[nodeListi]->numInternalElements();
      const auto rin = 2.0/vol[nodeListi]->nodeListPtr()->nodesPerSmoothingScale();

#pragma omp parallel for                                        \
  private(pairPlanes, voli)
      for (auto i = 0; i < n; ++i) {
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  rhoi = rho(nodeListi, i);
        auto        gradRhoi = gradRho(nodeListi, i);
        const auto  grhat = gradRhoi.unitVector();
        const auto  Hinvi = Hi.Inverse();
        const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        auto&       celli = polycells(nodeListi, i);
        const bool  interior = (surfacePoint(nodeListi, i) == 0);

        // First our own void points.
        pairPlanes.clear();
        for (const auto& etaVoid: etaVoidPoints(nodeListi, i)) {
          const auto rji = Hinvi*etaVoid;
          const auto nhat = -rji.unitVector();
          pairPlanes.push_back(Plane(rji, nhat));
        }

        // Add any neighbor void points this cell might interact with.
        for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          for (auto j: fullConnectivity[nodeListj]) {
            const auto& rj = position(nodeListj, j);
            const auto& Hj = H(nodeListj, j);
            const auto  Hinvj = Hj.Inverse();
            for (const auto& etaVoid: etaVoidPoints(nodeListj, j)) {
              const auto rji = Hinvj*etaVoid + 0.5*(rj - ri);
              const auto nhat = -rji.unitVector();
              pairPlanes.push_back(Plane(rji, nhat));
            }
          }
        }

        // Clip by the void neighbors.
        std::sort(pairPlanes.begin(), pairPlanes.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
        PolyClipper::clipPolygon(celli, pairPlanes);
        CHECK(not celli.empty());

        // Compute the moments of the clipped cell.
        PolyClipper::moments(voli, deltaMedian(nodeListi, i), celli);

        if (interior) {

          // We only use the volume result if interior.
          vol(nodeListi, i) = voli;

          // // Apply the gradient limiter;
          // gradRhoi *= phi;

          // Is there a significant density gradient?
          if (sqrt(gradRhoi.magnitude2()*vol(nodeListi, i)) >= 1e-8*rhoi) {

            const auto nhat1 = gradRhoi.unitVector();
            double x1, x2;
            findPolygonExtent(x1, x2, nhat1, celli);
            CHECK2(x1 <= 0.0 and x2 >= 0.0, nodeListi << " " << i << " " << ri << " " << x1 << " " << x2);
            const Scalar b = gradRhoi.magnitude();

            // This version uses the medial position.
            const auto thpt = sqrt(abs(rhoi*rhoi + rhoi*b*(x1 + x2) + 0.5*b*b*(x1*x1 + x2*x2)));
            const auto xm1 = -(rhoi + thpt)/b;
            const auto xm2 = (-rhoi + thpt)/b;
            const auto deltaCentroidi = deltaMedian(nodeListi, i);
            if (xm1 >= x1 and xm1 <= x2) {
              deltaMedian(nodeListi, i) = xm1*nhat1 - deltaCentroidi.dot(nhat1)*nhat1 + deltaCentroidi;
            } else {
              deltaMedian(nodeListi, i) = xm2*nhat1 - deltaCentroidi.dot(nhat1)*nhat1 + deltaCentroidi;
            }

            // // This version simply tries rho^2 weighting.
            // deltaMedian(nodeListi, i) = ((0.5*rhoi*(x2*x2 - x1*x1) +
            //                               2.0/3.0*rhoi*b*(x2*x2*x2 - x1*x1*x1) +
            //                               0.25*b*b*(x2*x2*x2*x2 - x1*x1*x1*x1))/
            //                              (pow3(rhoi + b*x2) - pow3(rhoi + b*x1)/(3.0*b)))*nhat1 - deltaCentroidi.dot(nhat1)*nhat1 + deltaCentroidi;
          }

          // OK, this is an interior point from the perspective that it was clipped within our critical
          // radius on all sides.  However, if we have a bounding polygon we may still want to call it a
          // surface if in fact there are still facets from that bounding polygon on this cell.
          // if (haveFacetedBoundaries and returnSurface) {
          //   auto vitr = celli.begin();
          //   while (interior and vitr != celli.end()) {
          //     interior = not pointOnPolygon(ri + vitr->position,
          //                                   factedBoundaries[nodeListi].vertices(),
          //                                   1.0e-8);
          //     ++vitr;
          //   }

          // }
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

        // Flip the surface bit if necesary.
        if (not interior) {
          // This is a point that touches the bounding polygon.  Flag it as surface.
          if (returnSurface) surfacePoint(nodeListi, i) |= 1;
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
