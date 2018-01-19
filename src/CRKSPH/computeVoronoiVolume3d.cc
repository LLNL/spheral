//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include <algorithm>
#include <utility>
#include <ctime>

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

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

using namespace FastMath;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NeighborSpace::Neighbor;
using NeighborSpace::ConnectivityMap;

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
                     const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& rho,
                     const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& gradRho,
                     const ConnectivityMap<Dim<3> >& connectivityMap,
                     const std::vector<Dim<3>::FacetedVolume>& boundaries,
                     const std::vector<std::vector<Dim<3>::FacetedVolume> >& holes,
                     const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                     const FieldList<Dim<3>, int>& voidPoint,
                     FieldList<Dim<3>, int>& surfacePoint,
                     FieldList<Dim<3>, Dim<3>::Scalar>& vol,
                     FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& deltaMedian,
                     FieldSpace::FieldList<Dim<3>, vector<Dim<3>::Vector>>& etaVoidPoints,
                     FieldSpace::FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells) {

  TIME_computeVoronoiVolume3d.start();

#ifdef NOR3D
  VERIFY2(false, "ERROR: computeVoronoiVolume requires compilation with R3D third party library.");
#else

  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef PolyClipper::Plane3d Plane;

  const auto numGens = position.numNodes();
  const auto numNodeLists = position.size();
  const auto numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());
  const auto numBounds = boundaries.size();
  const auto haveBoundaries = numBounds == numNodeLists;
  const auto haveWeights = weight.size() == numNodeLists;
  const auto returnSurface = surfacePoint.size() == numNodeLists;
  const auto returnCells = cells.size() == numNodeLists;

  REQUIRE(numBounds == 0 or numBounds == numNodeLists);
  REQUIRE(holes.size() == numBounds);

  // std::clock_t t0, 
  //   ttotal = std::clock_t(0), 
  //   tplanesneighbors = std::clock_t(0), 
  //   tplanesboundaries = std::clock_t(0), 
  //   tplanesort = std::clock_t(0), 
  //   tclip = std::clock_t(0), 
  //   tinterior = std::clock_t(0),
  //   tcentroid = std::clock_t(0),
  //   tsurface = std::clock_t(0),
  //   tbound = std::clock_t(0),
  //   tcell = std::clock_t(0);

  // ttotal = std::clock();

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

    // Walk the points.
    vector<Plane> pairPlanes, voidPlanes;
    double voli;
    PolyClipper::Polyhedron celli;
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = vol[nodeListi]->numInternalElements();
      const auto rin = 2.0/vol[nodeListi]->nodeListPtr()->nodesPerSmoothingScale();
#pragma omp parallel for                        \
  private(pairPlanes, voidPlanes, voli, celli)
      for (auto i = 0; i < n; ++i) {

        const auto& ri = position(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  rhoi = rho(nodeListi, i);
        auto        gradRhoi = gradRho(nodeListi, i);
        const auto  grhat = gradRhoi.unitVector();
        const auto  Hdeti = Hi.Determinant();
        const auto  Hinv = Hi.Inverse();
        const auto  weighti = haveWeights ? weight(nodeListi, i) : 1.0;
        const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

        // Prepare to accumulate any void point positions.
        pairPlanes.clear();
        voidPlanes.clear();

        // t0 = std::clock();

        // Initialize our seed cell shape.
        celli = cell0;
        for (auto& v: celli) v.position = 1.1*rin*Hinv*v.position;

        // Clip by any boundaries first.
        if (haveBoundaries) {
          const auto& facets = boundaries[nodeListi].facets();
          CHECK(boundaries[nodeListi].contains(ri, false));
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

        // Add clipping planes from neighbors.
        for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

          // Check if this point has multi-material neighbors.
          if (returnSurface and 
              nodeListj != nodeListi and 
              not fullConnectivity[nodeListj].empty()) surfacePoint(nodeListi, i) |= (1 << (nodeListj + 1));

          for (auto jItr = fullConnectivity[nodeListj].begin();
               jItr != fullConnectivity[nodeListj].end();
               ++jItr) {
            const auto  j = *jItr;
            const auto& rj = position(nodeListj, j);
            const auto  weightj = haveWeights ? weight(nodeListj, j) : 1.0;

            // Build the planes for the clipping half-spaces.
            const auto rji = rj - ri;
            const auto nhat = -rji.unitVector();
            const auto wij = weighti/(weightj + weighti);
            if (voidPoint(nodeListj, j) == 0) {
              pairPlanes.push_back(Plane(wij*rji, nhat));
            } else {
              voidPlanes.push_back(Plane(wij*rji, nhat));
            }
          }
        }

        // tplanesneighbors += std::clock() - t0;
        // t0 = std::clock();

        // Sort the planes by distance -- let's us clip more efficiently.
        std::sort(pairPlanes.begin(), pairPlanes.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });

        // tplanesort += std::clock() - t0;
        // t0 = std::clock();

        // Clip by non-void neighbors first.
        PolyClipper::clipPolyhedron(celli, pairPlanes);
        CHECK(celli.size() > 0);

        // tclip += std::clock() - t0;

        // Check if the final polyghedron is entirely within our "interior" check radius.
        // We preserve remaining original vertices as void points.
        bool interior = true;
        // t0 = std::clock();
        {
          for (const auto& vert: celli) {
            const auto peta = Hi*vert.position;
            if (peta.magnitude2() > rin*rin) {
              interior = false;
              if (returnSurface) {
                surfacePoint(nodeListi, i) |= 1;
                const Vector etaj = 0.5*rin*peta.unitVector();
                etaVoidPoints(nodeListi, i).push_back(etaj);
                const auto rji = ri - Hinv*etaj;
                const auto nhat = -rji.unitVector();
                voidPlanes.push_back(Plane(0.5*rji, nhat));
              }
            }
          }
        }
        // tinterior += std::clock() - t0;

        // Clip the cell geometry by the void planes.
        std::sort(voidPlanes.begin(), voidPlanes.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
        PolyClipper::clipPolyhedron(celli, voidPlanes);
        CHECK(celli.size() > 0);

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
          // tcentroid += std::clock() - t0;

          // OK, this is an interior point from the perspective that it was clipped within our critical
          // radius on all sides.  However, if we have a bounding polygon we may still want to call it a
          // surface if in fact there are still facets from that bounding polygon on this cell.
          // t0 = std::clock();
          if (haveBoundaries and returnSurface) {
            auto vitr = celli.begin();
            while (interior and vitr != celli.end()) {
              interior = not pointOnPolyhedron(ri + vitr->position,
                                               boundaries[nodeListi].vertices(),
                                               1.0e-8);
              ++vitr;
            }
          }
          // tsurface += std::clock() - t0;
        }

        // Check if the candidate motion is still in the boundary.  If not, project back.
        // t0 = std::clock();
        if (haveBoundaries) {
          if (not boundaries[nodeListi].contains(ri + deltaMedian(nodeListi, i), false)) {
            deltaMedian(nodeListi, i) = boundaries[nodeListi].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
          }
          for (unsigned ihole = 0; ihole != holes[nodeListi].size(); ++ihole) {
            if (holes[nodeListi][ihole].contains(ri + deltaMedian(nodeListi, i), false)) {
              deltaMedian(nodeListi, i) = holes[nodeListi][ihole].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
            }
          }
        }
        // tbound += std::clock() - t0;

        // Flip the surface bit if necesary.
        // if (not interior) {
        //   // This is a point that touches the bounding polygon.  Flag it as surface.
        //   if (returnSurface) surfacePoint(nodeListi, i) |= 1;
        // }
        // tbound += std::clock() - t0;

        // If requested, we can return the cell geometries.
        if (returnCells) {
          // t0 = std::clock();
          // PolyClipper::collapseDegenerates(celli, 1.0e-10);
          PolyClipper::convertFromPolyhedron(cells(nodeListi, i), celli);
          cells(nodeListi, i) += ri;
          // tcell += std::clock() - t0;
        }
      }
    }
  }

  // ttotal = std::clock() - ttotal;
  // if (Process::getRank() == 0) cout << "computeVoronoiVolume3d timing: " 
  //                                   << "tplanesneighbors=" << (tplanesneighbors / (double) CLOCKS_PER_SEC) 
  //                                   << " tplanesboundaries=" << (tplanesboundaries / (double) CLOCKS_PER_SEC) 
  //                                   << " tplanesort=" << (tplanesort / (double) CLOCKS_PER_SEC) 
  //                                   << " tclip=" << (tclip / (double) CLOCKS_PER_SEC) 
  //                                   << " tinterior=" << (tinterior / (double) CLOCKS_PER_SEC) 
  //                                   << " tcentroid=" << (tcentroid / (double) CLOCKS_PER_SEC) 
  //                                   << " tsurface=" << (tsurface / (double) CLOCKS_PER_SEC) 
  //                                   << " tbound=" << (tbound / (double) CLOCKS_PER_SEC) 
  //                                   << " tcell=" << (tcell / (double) CLOCKS_PER_SEC) 
  //                                   << " ttotal=" << (ttotal / (double) CLOCKS_PER_SEC) << endl;
#endif

  TIME_computeVoronoiVolume3d.stop();
}
    
}
}
