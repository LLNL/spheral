//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/pointOnPolyhedron.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/Timer.hh"
#include "Geometry/polyclipper.hh"

extern Timer TIME_computeVoronoiVolume3d;

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

//------------------------------------------------------------------------------
// Find the 1D extent of a polygon along the given direction.
//------------------------------------------------------------------------------
inline
void findPolyhedronExtent(double& xmin, double& xmax, 
                          const Dim<3>::Vector& nhat, 
                          const PolyClipper::Polyhedron& celli) {
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
// 3D
//------------------------------------------------------------------------------
void
computeVoronoiVolume(const FieldList<Dim<3>, Dim<3>::Vector>& position,
                     const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                     const FieldList<Dim<3>, Dim<3>::Scalar>& rho,
                     const FieldList<Dim<3>, Dim<3>::Vector>& gradRho,
                     const ConnectivityMap<Dim<3> >& connectivityMap,
                     const FieldList<Dim<3>, Dim<3>::SymTensor>& damage,
                     const std::vector<Dim<3>::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<Dim<3>::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dim<3>>*>& boundaries,
                     const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                     const FieldList<Dim<3>, int>& voidPoint,
                     FieldList<Dim<3>, int>& surfacePoint,
                     FieldList<Dim<3>, Dim<3>::Scalar>& vol,
                     FieldList<Dim<3>, Dim<3>::Vector>& deltaMedian,
                     FieldList<Dim<3>, vector<Dim<3>::Vector>>& etaVoidPoints,
                     FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells) {

  TIME_computeVoronoiVolume3d.start();

  // Pre-conditions
  REQUIRE(facetedBoundaries.size() == 0 or facetedBoundaries.size() == position.size());
  REQUIRE(holes.size() == facetedBoundaries.size());

  typedef Dim<3> Dimension;
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef PolyClipper::Plane3d Plane;

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

    // Build an approximation of the starting kernel shape (in eta space) as an icosahedron with vertices
    // at a unit radius.
    PolyClipper::Polyhedron cell0;
    {
      const auto t = (1.0 + sqrt(5.0)) / 2.0;
      const vector<Vector> vertsIco = {           // Array of vertex coordinates.
        Vector(-1,  t,  0)/sqrt(1 + t*t),
        Vector( 1,  t,  0)/sqrt(1 + t*t),
        Vector(-1, -t,  0)/sqrt(1 + t*t),
        Vector( 1, -t,  0)/sqrt(1 + t*t),
        Vector( 0, -1,  t)/sqrt(1 + t*t),
        Vector( 0,  1,  t)/sqrt(1 + t*t),
        Vector( 0, -1, -t)/sqrt(1 + t*t),
        Vector( 0,  1, -t)/sqrt(1 + t*t),
        Vector( t,  0, -1)/sqrt(1 + t*t),
        Vector( t,  0,  1)/sqrt(1 + t*t),
        Vector(-t,  0, -1)/sqrt(1 + t*t),
        Vector(-t,  0,  1)/sqrt(1 + t*t)
      };
      const vector<vector<unsigned>> facesIco = {
        // 5 faces around point 0
        {0, 11, 5},
        {0, 5, 1},
        {0, 1, 7},
        {0, 7, 10},
        {0, 10, 11},
        // 5 adjacent faces
        {1, 5, 9},
        {5, 11, 4},
        {11, 10, 2},
        {10, 7, 6},
        {7, 1, 8},
        // 5 faces around point 3
        {3, 9, 4},
        {3, 4, 2},
        {3, 2, 6},
        {3, 6, 8},
        {3, 8, 9},
        // 5 adjacent faces
        {4, 9, 5},
        {2, 4, 11},
        {6, 2, 10},
        {8, 6, 7},
        {9, 8, 1}
      };
      PolyClipper::convertToPolyhedron(cell0, FacetedVolume(vertsIco, facesIco));
    }
    const auto nverts = cell0.size();

    // We'll need to hang onto the PolyClipper cells.
    FieldList<Dim<3>, PolyClipper::Polyhedron> polycells(FieldStorageType::CopyFields);
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      polycells.appendNewField("polycells", vol[nodeListi]->nodeList(), PolyClipper::Polyhedron());
    }

    // First pass: clip by neighbors and generate void points.
    double vol0, voli;
    vector<Plane> pairPlanes;
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = vol[nodeListi]->numInternalElements();
      const auto rin = 2.0/vol[nodeListi]->nodeListPtr()->nodesPerSmoothingScale();

#pragma omp parallel for                \
  private(pairPlanes, vol0, voli)
      for (auto i = 0; i < n; ++i) {
        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hinv = Hi.Inverse();
        const auto  weighti = haveWeights ? weight(nodeListi, i) : 1.0;
        const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        auto&       celli = polycells(nodeListi, i);

        // Initialize our seed cell shape.
        celli = cell0;
        for (auto& v: celli) v.position = 1.1*rin*Hinv*v.position;

        // Clip by any faceted boundaries first.
        pairPlanes.clear();
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
        PolyClipper::clipPolyhedron(celli, pairPlanes);
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
            PolyClipper::clipPolyhedron(celli, pairPlanes);
            if (returnSurface) {
              vol0 = voli;
              PolyClipper::moments(voli, deltaMedian(nodeListi, i), celli);
              if (voli < vol0) surfacePoint(nodeListi, i) |= (1 << (nodeListj + 1));
            }
          }
        }

        // Check if the final polyghedron is entirely within our "interior" check radius.
        // We preserve remaining original vertices as void points.
        for (const auto& vert: celli) {
          const auto peta = Hi*vert.position;
          if (peta.magnitude2() > rin*rin) {
            if (returnSurface) {
              surfacePoint(nodeListi, i) |= 1;
              const Vector etaj = 0.5*rin*peta.unitVector();
              etaVoidPoints(nodeListi, i).push_back(etaj);
            }
          }
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
        bool        interior = (surfacePoint(nodeListi, i) == 0);

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
        PolyClipper::clipPolyhedron(celli, pairPlanes);
        CHECK(not celli.empty());

        // Compute the moments of the clipped cell.
        PolyClipper::moments(voli, deltaMedian(nodeListi, i), celli);

        if (interior) {

          // We only use the volume result if interior.
          vol(nodeListi, i) = voli;

          // // Apply the gradient limiter;
          // gradRhoi *= phi;

          // Is there a significant density gradient?
          if (gradRhoi.magnitude()*Dim<3>::rootnu(vol(nodeListi, i)) >= 1e-8*rhoi) {

            const auto nhat1 = gradRhoi.unitVector();
            double x1, x2;
            findPolyhedronExtent(x1, x2, nhat1, celli);
            CHECK2(x1 <= 0.0 and x2 >= 0.0, nodeListi << " " << i << " " << ri << " " << x1 << " " << x2);
            const Scalar b = gradRhoi.magnitude();

            // This version uses the medial position.
            const Scalar thpt = sqrt(abs(rhoi*rhoi + rhoi*b*(x1 + x2) + 0.5*b*b*(x1*x1 + x2*x2)));
            const Scalar xm1 = -(rhoi + thpt)/b;
            const Scalar xm2 = (-rhoi + thpt)/b;
            const auto deltaCentroidi = deltaMedian(nodeListi, i);
            if (xm1 >= x1 and xm1 <= x2) {
              deltaMedian(nodeListi, i) = xm1*nhat1 - deltaCentroidi.dot(nhat1)*nhat1 + deltaCentroidi;
            } else {
              deltaMedian(nodeListi, i) = xm2*nhat1 - deltaCentroidi.dot(nhat1)*nhat1 + deltaCentroidi;
            }
          }

          // OK, this is an interior point from the perspective that it was clipped within our critical
          // radius on all sides.  However, if we have a bounding polygon we may still want to call it a
          // surface if in fact there are still facets from that bounding polygon on this cell.
          if (haveFacetedBoundaries and returnSurface) {
            auto vitr = celli.begin();
            while (interior and vitr != celli.end()) {
              interior = not pointOnPolyhedron(ri + vitr->position,
                                               facetedBoundaries[nodeListi],
                                               1.0e-8);
              ++vitr;
            }
          }
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
        // if (not interior) {
        //   // This is a point that touches the bounding polygon.  Flag it as surface.
        //   if (returnSurface) surfacePoint(nodeListi, i) |= 1;
        // }

        // If requested, we can return the cell geometries.
        if (returnCells) {
          // PolyClipper::collapseDegenerates(celli, 1.0e-5);
          PolyClipper::convertFromPolyhedron(cells(nodeListi, i), celli);
          cells(nodeListi, i) += ri;
        }
      }
    }
  }

  TIME_computeVoronoiVolume3d.stop();
}
    
}
