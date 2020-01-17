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

extern Timer TIME_computeVoronoiVolume;

#include <algorithm>
#include <utility>
#include <limits>
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

namespace {  // anonymous namespace
//------------------------------------------------------------------------------
// Trait class for local Dimension
//------------------------------------------------------------------------------
template<typename Dimension> struct ClippingType;

//..............................................................................
// 2D
template<>
struct
ClippingType<Dim<2>> {
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef PolyClipper::Plane2d Plane;
  typedef PolyClipper::Polygon PolyVolume;
  
  // Unit circle as template shape.
  static PolyVolume unitPolyVolume() {
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
    return cell0;
  }

  // clipping operation
  static void clip(PolyVolume& cell, const std::vector<Plane>& planes) {
    PolyClipper::clipPolygon(cell, planes);
  }

  // moment computation
  static void moments(double& vol, Vector& cent, const PolyVolume& cell) {
    PolyClipper::moments(vol, cent, cell);
  }

  // collapse degenerate points
  static void collapseDegenerates(PolyVolume& cell, const double tol) {
    PolyClipper::collapseDegenerates(cell, tol);
  }

  // Convert PolyClipper::Polygon -> Spheral::Polygon
  static std::vector<std::set<int>> convertFromPolyVolume(FacetedVolume& spheralcell, const PolyVolume& polycell) {
    return PolyClipper::convertFromPolygon(spheralcell, polycell);
  }

  // Generate the reduced void point stencil -- up to 4 for 2D
  static std::vector<Vector> createEtaVoidPoints(const Vector& etaVoidAvg,
                                                 const int nvoid,
                                                 const double rin,
                                                 const SymTensor& Hi,
                                                 const SymTensor& Hinvi,
                                                 const PolyVolume& celli) {
    std::vector<Vector> result;
    const auto nverts = 18;
    const auto thetaVoidAvg = atan2(etaVoidAvg.y(), etaVoidAvg.x());
    const auto nv = max(1U, min(4U, unsigned(4.0*double(nvoid)/double(nverts))));
    for (unsigned k = 0; k != nv; ++k) {
      const auto theta = thetaVoidAvg + (0.5*k - 0.25*(nv - 1))*M_PI;
      const auto etaVoid = Vector(0.5*rin*cos(theta), 0.5*rin*sin(theta));
      result.push_back(etaVoid);
      const auto rji = Hinvi*etaVoid;
      const auto nhat = -rji.unitVector();
    }
    ENSURE(result.size() == nv);
    return result;
  }
    
  // toString
  static std::string toString(const PolyVolume& celli) {
    return PolyClipper::polygon2string(celli);
  }

};

//..............................................................................
// 3D
template<>
struct
ClippingType<Dim<3>> {
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef PolyClipper::Plane3d Plane;
  typedef PolyClipper::Polyhedron PolyVolume;
  
  // Build an approximation of the starting kernel shape (in eta space) as an icosahedron with vertices
  static PolyVolume unitPolyVolume() {
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
    PolyClipper::Polyhedron cell0;
    PolyClipper::convertToPolyhedron(cell0, FacetedVolume(vertsIco, facesIco));
    ENSURE(cell0.size() == 12);
    return cell0;
  }

  // clipping operation
  static void clip(PolyVolume& cell, const std::vector<Plane>& planes) {
    PolyClipper::clipPolyhedron(cell, planes);
  }

  // moment computation
  static void moments(double& vol, Vector& cent, const PolyVolume& cell) {
    PolyClipper::moments(vol, cent, cell);
  }

  // collapse degenerate points
  static void collapseDegenerates(PolyVolume& cell, const double tol) {
    PolyClipper::collapseDegenerates(cell, tol);
  }

  // Convert PolyClipper::Polyhedron -> Spheral::Polyhedron
  static std::vector<std::set<int>> convertFromPolyVolume(FacetedVolume& spheralcell, const PolyVolume& polycell) {
    return PolyClipper::convertFromPolyhedron(spheralcell, polycell);
  }

  // In 3D we simply use any unclipped original vertices as void generators
  static std::vector<Vector> createEtaVoidPoints(const Vector& etaVoidAvg,
                                                 const int nvoid,
                                                 const double rin,
                                                 const SymTensor& Hi,
                                                 const SymTensor& Hinvi,
                                                 const PolyVolume& celli) {
    std::vector<Vector> result;
    for (const auto& vert: celli) {
      const auto peta = Hi*vert.position;
      if (peta.magnitude2() > rin*rin) {
        const Vector etaj = 0.5*rin*peta.unitVector();
        result.push_back(etaj);
      }
    }
    return result;
  }
    
  // toString
  static std::string toString(const PolyVolume& celli) {
    return PolyClipper::polyhedron2string(celli);
  }

};

//------------------------------------------------------------------------------
// Return the set of face flags for surface points
// We come in with the following convention for the clips on the vertices:
//   1.  Clipped by a bounding Poly* (exterior or interior): numeric_limits<int>::min()
//   2.  Clipped by a Node in another NodeList             : ~nodeListj
//   3.  Clipped by a void plane                           : >=0
//------------------------------------------------------------------------------
// 2D
std::vector<CellFaceFlag> extractFaceFlags(const GeomPolygon& cell,
                                           const std::vector<std::set<int>>& vertexClips) {
  REQUIRE(vertexClips.size() == cell.vertices().size());
  const auto boundingSurfaceClipFlag = std::numeric_limits<int>::min();    // Special cellFaceFlag to indicate clipping by a bounding polyhedron/hole.
  const auto& facets = cell.facets();
  const auto  nfacets = facets.size();
  std::vector<CellFaceFlag> result;
  for (auto ifacet = 0; ifacet < nfacets; ++ifacet) {
    const auto& facet = facets[ifacet];
    const auto& clips1 = vertexClips[facet.ipoint1()];
    const auto& clips2 = vertexClips[facet.ipoint2()];
    // if (not (clips1.empty() and clips2.empty())) {
    //   cerr << " --> " << ifacet << " " << facet.ipoint1() << " " << facet.ipoint2() << endl
    //        << "    ";
    //   for (auto x: clips1) cerr << " " << x;
    //   cerr << endl
    //        << "    ";
    //   for (auto x: clips2) cerr << " " << x;
    //   cerr << endl;
    // }
    if (clips1.size() > 0 and clips2.size() > 0) {
      for (auto iplane1: clips1) {
        if (clips2.find(iplane1) != clips2.end()) {
          if (iplane1 == boundingSurfaceClipFlag or
              iplane1 >= 0) {
            result.push_back(CellFaceFlag({ifacet, -1, -1}));       // Clipped by a boundary/void
          } else {
            result.push_back(CellFaceFlag({ifacet, ~iplane1, -1})); // Clipped by another NodeList
          }
        }
        break;
      }
    }
  }
  return result;
}

// 3D
std::vector<CellFaceFlag> extractFaceFlags(const GeomPolyhedron& cell,
                                           const std::vector<std::set<int>>& vertexClips) {
  REQUIRE(vertexClips.size() == cell.vertices().size());
  const auto boundingSurfaceClipFlag = std::numeric_limits<int>::min();    // Special cellFaceFlag to indicate clipping by a bounding polyhedron/hole.
  const auto& facets = cell.facets();
  const auto  nfacets = facets.size();
  std::vector<CellFaceFlag> result;
  for (auto ifacet = 0; ifacet < nfacets; ++ifacet) {
    const auto& facet = facets[ifacet];
    const auto& ipoints = facet.ipoints();
    const auto  npoints = ipoints.size();
    auto intersection = vertexClips[ipoints[0]];
    auto k = 1;
    while (not intersection.empty() and k < npoints) {
      std::set<int> newIntersection;
      std::set_intersection(intersection.begin(), intersection.end(),
                            vertexClips[ipoints[k]].begin(), vertexClips[ipoints[k]].end(),
                            std::inserter(newIntersection, newIntersection.begin()));
      intersection = newIntersection;
      ++k;
    }
    if (not intersection.empty()) {
      CHECK(intersection.size() == 1);
      const auto iplane1 = *intersection.begin();
      if (iplane1 == boundingSurfaceClipFlag or
          iplane1 >= 0) {
        result.push_back(CellFaceFlag({ifacet, -1, -1}));       // Clipped by a boundary/void
      } else {
        result.push_back(CellFaceFlag({ifacet, ~iplane1, -1})); // Clipped by another NodeList
      }
    }
  }
  return result;
}

}           // anonymous namespace

//------------------------------------------------------------------------------
// The method in question
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeVoronoiVolume(const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     const ConnectivityMap<Dimension >& connectivityMap,
                     const FieldList<Dimension, typename Dimension::SymTensor>& damage,
                     const std::vector<typename Dimension::FacetedVolume>& facetedBoundaries,
                     const std::vector<std::vector<typename Dimension::FacetedVolume> >& holes,
                     const std::vector<Boundary<Dimension>*>& boundaries,
                     const FieldList<Dimension, typename Dimension::Scalar>& weight,
                     FieldList<Dimension, int>& surfacePoint,
                     FieldList<Dimension, typename Dimension::Scalar>& vol,
                     FieldList<Dimension, typename Dimension::Vector>& deltaMedian,
                     FieldList<Dimension, vector<typename Dimension::Vector>>& etaVoidPoints,
                     FieldList<Dimension, typename Dimension::FacetedVolume>& cells,
                     FieldList<Dimension, std::vector<CellFaceFlag>>& cellFaceFlags) {

  TIME_computeVoronoiVolume.start();

  // Pre-conditions
  REQUIRE(facetedBoundaries.size() == 0 or facetedBoundaries.size() == position.size());
  REQUIRE(vol.size() == position.size());
  REQUIRE(deltaMedian.size() == position.size());
  REQUIRE(holes.size() == facetedBoundaries.size());

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef typename ClippingType<Dimension>::Plane Plane;
  typedef typename ClippingType<Dimension>::PolyVolume PolyVolume;

  const auto numGens = position.numNodes();
  const auto numNodeLists = position.size();
  const auto numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());
  const auto haveFacetedBoundaries = facetedBoundaries.size() == numNodeLists;
  const auto haveBoundaries = not boundaries.empty();
  const auto haveWeights = weight.size() == numNodeLists;
  const auto haveDamage = false;  // damage.size() == numNodeLists;   // Suspending the idea of forcing surface based on damage
  const auto returnSurface = surfacePoint.size() == numNodeLists;
  const auto returnCells = cells.size() == numNodeLists;
  const auto returnCellFaceFlags = cellFaceFlags.size() == numNodeLists;

  const auto boundingSurfaceClipFlag = std::numeric_limits<int>::min();    // Special cellFaceFlag to indicate clipping by a bounding polyhedron/hole.

  if (returnSurface) {
    surfacePoint = 0;
    etaVoidPoints = vector<Vector>();
  }

  if (numGensGlobal > 0) {

    const auto& pairs = connectivityMap.nodePairList();
    const auto  npairs = pairs.size();

    // The criterion for the initial starting sphere of a point in eta space.
    const auto rin = 2.0/vol[0]->nodeListPtr()->nodesPerSmoothingScale();
    // const auto etaMax = vol[0]->nodeList().neighbor().kernelExtent();
    // const auto rin = etaMax;

    // Unit circle as template shape.
    auto cell0 = ClippingType<Dimension>::unitPolyVolume();

    // We'll need to hang onto the PolyClipper cells and any per cell void points.
    FieldList<Dimension, PolyVolume> polycells(FieldStorageType::CopyFields);
    FieldList<Dimension, vector<vector<Plane>>> pairPlanes(FieldStorageType::CopyFields);
    FieldList<Dimension, vector<Plane>> voidPlanes(FieldStorageType::CopyFields);
    for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      polycells.appendNewField("polycells", vol[nodeListi]->nodeList(), cell0);
      pairPlanes.appendNewField("pair planes", vol[nodeListi]->nodeList(), vector<vector<Plane>>(numNodeLists));
      voidPlanes.appendNewField("void planes", vol[nodeListi]->nodeList(), vector<Plane>());
    }

#pragma omp parallel
    {

      //==========================================================================
      // First pass: clip by any faceted boundaries/holes.
      // cerr << "FIRST pass after polyhedral boundary clipping" << endl;
      for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
        const auto ni = polycells[nodeListi]->numInternalElements();
#pragma omp for
        for (auto i = 0; i < ni; ++i) {

          // Initialize the per point polygon by scaling by its H extent
          const auto& ri = position(nodeListi, i);
          const auto& Hi = H(nodeListi, i);
          const auto  Hinv = Hi.Inverse();
#pragma omp critical (computeVoronoiVolume_polycells)
          {
            for (auto& v: polycells(nodeListi, i)) v.position = 1.1*rin*Hinv*v.position;
          }

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
#pragma omp critical (computeVoronoiVolume_polycells)
            {
              ClippingType<Dimension>::clip(polycells(nodeListi, i), boundPlanes);

              // Did we make any new (and therefore external) surfaces on the cell?
              // If so, mark those vertices as external.
              for (auto& vertex: polycells(nodeListi, i)) {
                if (not vertex.clips.empty()) vertex.clips = {boundingSurfaceClipFlag};
              }
            }
          }
        }
      }

      //==========================================================================
      // Second pass: clip by neighbor points.  Note we have to keep track of
      // which NodeLists actually clip each polygon in order to detect material
      // surfaces.
      // cerr << "SECOND pass after node-node clipping" << endl;

      // Thread private scratch variables
      int i, j, nodeListi, nodeListj;
      auto pairPlanes_thread = pairPlanes.threadCopy(ThreadReduction::SUM, true);  // force copying the original FieldList

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
        CHECK(pairPlanesi.size() == numNodeLists);

        // State of node j
        const auto& rj = position(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  Hjnv = Hj.Inverse();
        const auto  weightj = haveWeights ? weight(nodeListj, j) : 1.0;
        auto&       pairPlanesj = pairPlanes_thread(nodeListj, j);
        CHECK(pairPlanesj.size() == numNodeLists);

        // Build the planes for the clipping half-spaces.
        const auto rji = rj - ri;
        const auto nhat = -rji.unitVector();
        const auto wij = weighti/(weighti + weightj);
        const auto wji = weightj/(weighti + weightj);
        pairPlanesi[nodeListj].push_back(Plane( wij*rji,  nhat, kk));
        pairPlanesj[nodeListi].push_back(Plane(-wji*rji, -nhat, kk));
      }

      // Collect the pair planes across threads.
#pragma omp critical (computeVoronoiVolume_pass2)
      {
        pairPlanes_thread.threadReduce();
      }
// #pragma omp barrier  // Wait for the pairPlanes to be done
    }   // OMP parallel

#pragma omp parallel
    {
      // Clip by the neighbors, and look for any locally generated void points.
      auto etaVoidPoints_thread = etaVoidPoints.threadCopy(); // ThreadReduction::SUM, true);
      PolyVolume celli;
      for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
        const auto ni = etaVoidPoints_thread[nodeListi]->numInternalElements();
#pragma omp parallel for
        for (auto i = 0; i < ni; ++i) {
          const auto& Hi = H(nodeListi, i);
          const auto  Hinvi = Hi.Inverse();
          auto        pairPlanesi = pairPlanes(nodeListi, i);  // Deliberately make a copy
#pragma omp critical (computeVornoiVolume_polycells)
          {
            celli = polycells(nodeListi, i);         // Also make a copy the starting global cell for this point
          }
          CHECK(not celli.empty());

          if (returnSurface) {
            // If we're returning the surface we have to pay attention to which materials actually clip this point.
            double vol0, vol1;
            Vector cent;
            ClippingType<Dimension>::moments(vol0, cent, celli);
            for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
              std::sort(pairPlanesi[nodeListj].begin(), pairPlanesi[nodeListj].end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
              ClippingType<Dimension>::clip(celli, pairPlanesi[nodeListj]);
              ClippingType<Dimension>::moments(vol1, cent, celli);
              if (vol1 < vol0) {
                vol0 = vol1;
                if (nodeListj != nodeListi) {
                  surfacePoint(nodeListi, i) |= (1 << (nodeListj + 1));
                  for (auto& vertex: celli) {
                    std::set<int> newclips;
                    for (const auto iplane: vertex.clips) {
                      if (iplane == boundingSurfaceClipFlag) {
                        newclips.insert(boundingSurfaceClipFlag);
                      } else {
                        CHECK(iplane < pairPlanesi[nodeListj].size());
                        newclips.insert(~nodeListj);
                      }
                    }
                    vertex.clips = newclips;
                  }
                }
              }
            }
          } else {
            // Otherwise just clip by all the pair planes at once.
            for (auto nodeListj = 1; nodeListj < numNodeLists; ++nodeListj) pairPlanesi[0] += pairPlanesi[nodeListj];
            std::sort(pairPlanesi[0].begin(), pairPlanesi[0].end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
            ClippingType<Dimension>::clip(celli, pairPlanesi[0]);
          }
          CHECK(not celli.empty());

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
          auto& etaVoidPointsi = etaVoidPoints_thread(nodeListi, i);
          CHECK2(etaVoidPointsi.empty(), "(" << nodeListi << " " << i << ") : " << etaVoidPointsi.size());
          // Reduce the number of void points for this point to a reasonable stencil -- up to 4 for 2D.
          if (nvoid > 0) {
            etaVoidPointsi = ClippingType<Dimension>::createEtaVoidPoints(etaVoidAvg, nvoid, rin, Hi, Hinvi, celli);
          }

          // If this point is sufficiently damaged, we also create void points along the damaged directions.
          if (haveDamage and damage(nodeListi, i).Trace() > 1.0 - 1.0e-5) {
            const auto ev = damage(nodeListi, i).eigenVectors();
            for (auto jdim = 0; jdim < Dimension::nDim; ++jdim) {
              if (ev.eigenValues(jdim) > 1.0 - 1.0e-5) {
                const auto evecj = ev.eigenVectors.getColumn(jdim);
                etaVoidPointsi.push_back(-0.5*rin*evecj);
                etaVoidPointsi.push_back( 0.5*rin*evecj);
              }
            }
          }

          // Store the clipped cell thus far
#pragma omp critical (computeVornoiVolume_polycells)
          {
            polycells(nodeListi, i) = celli;
          }

        }   // end over i OMP parallel for
      }     // end over NodeLists

#pragma omp critical (computeVoronoiVolume_pass2_reduce)
      {
        etaVoidPoints_thread.threadReduce();
      }     // OMP critical
    }       // OMP parallel

    // Apply boundary conditions to the void points.
    if (not boundaries.empty()) {
      cerr << "Starting boundaries..." << endl;
      MPI_Barrier(Communicator::communicator());
      for (const auto& bc: boundaries) bc->applyFieldListGhostBoundary(etaVoidPoints);
      for (const auto& bc: boundaries) bc->finalizeGhostBoundary();
      MPI_Barrier(Communicator::communicator());
      cerr << "Finished boundaries..." << endl;
    }

#pragma omp parallel
    {
      //==========================================================================
      // Third pass: clip by any void points.
      // Start by adding any void clip planes from neighbors.
      // cerr << " --> " << omp_get_thread_num() << " THIRD PASS -- void clipping" << endl;
      int i, j, nodeListi, nodeListj;
      auto voidPlanes_thread = voidPlanes.threadCopy();
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
        const auto& etaVoidPointsi = etaVoidPoints(nodeListi, i);
        auto&       voidPlanesi = voidPlanes_thread(nodeListi, i);

        // State of node j
        const auto& rj = position(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  Hinvj = Hj.Inverse();
        const auto& etaVoidPointsj = etaVoidPoints(nodeListj, j);
        auto&       voidPlanesj = voidPlanes_thread(nodeListj, j);

        for (const auto etaVoid: etaVoidPoints(nodeListj, j)) {
          const auto rji = Hinvj*etaVoid + 0.5*(rj - ri);
          const auto nhat = -rji.unitVector();
          voidPlanesi.push_back(Plane(rji, nhat, kk));
        }

        for (const auto etaVoid: etaVoidPoints(nodeListi, i)) {
          const auto rij = Hinvi*etaVoid + 0.5*(ri - rj);
          const auto nhat = -rij.unitVector();
          voidPlanesj.push_back(Plane(rij, nhat, kk));
        }
      }

#pragma omp critical (computeVoronoiVolume_pass3_reduce_voidPlanes)
      {
        voidPlanes_thread.threadReduce();
      }
      // cerr << " --> " << omp_get_thread_num() << " THIRD PASS -- finished building voidPlanes" << endl;
#pragma omp barrier

      // Now we can do the void clipping, compute the final volumes, and finalize
      // surface detection.
      for (auto nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
        const auto n = vol[nodeListi]->numInternalElements();
#pragma omp parallel for
        for (auto i = 0; i < n; ++i) {
          const auto& ri = position(nodeListi, i);
          const auto& Hi = H(nodeListi, i);
          const auto  Hinvi = Hi.Inverse();
          auto        voidPlanesi = voidPlanes(nodeListi, i);  // Deliberate copy for thread safety
          PolyVolume celli;
#pragma omp critical (computeVoronoiVolume_polycells)
          {
            celli = polycells(nodeListi, i);         // Deliberate copy for thread safety
          }
          CHECK(not celli.empty());

          // Clip by the void planes, compute the volume, and check if the void surface had any effect.
          double vol0, vol1;
          std::sort(voidPlanesi.begin(), voidPlanesi.end(), [](const Plane& lhs, const Plane& rhs) { return lhs.dist < rhs.dist; });
          ClippingType<Dimension>::moments(vol0, deltaMedian(nodeListi, i), celli);
          ClippingType<Dimension>::clip(celli, voidPlanesi);
          CHECK(not celli.empty());
          ClippingType<Dimension>::moments(vol1, deltaMedian(nodeListi, i), celli);

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
            ClippingType<Dimension>::collapseDegenerates(celli, 1.0e-10);
#pragma omp critical (computeVoronoiVolume_pass3)
            {
              auto vertexClips = ClippingType<Dimension>::convertFromPolyVolume(cells(nodeListi, i), celli);
              cells(nodeListi, i) += ri;

              // If we're returning CellFaceFlags, build them.  Note -- we currently do not store which neighbor node
              // is responsible for each clipped facet.  Just the material or void.
              if (returnCellFaceFlags) {
                if (surfacePoint(nodeListi, i) == 0) {
                  cellFaceFlags(nodeListi, i).clear();
                } else {
                  cellFaceFlags(nodeListi, i) = extractFaceFlags(cells(nodeListi, i), vertexClips);
                }
              }
            }
          }
        }
      }
    }  // OMP parallel
  }    // numGensGlobal > 0

  cerr << "computeVoronoiVolume FINISHED" << endl;
  TIME_computeVoronoiVolume.stop();
}
    
//------------------------------------------------------------------------------
// Instantiations
//------------------------------------------------------------------------------
#ifdef SPHERAL2D
template
void
computeVoronoiVolume<Dim<2>>(const FieldList<Dim<2>, Dim<2>::Vector>& position,
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
                             FieldList<Dim<2>, std::vector<CellFaceFlag>>& cellFaceFlags);
#endif

#ifdef SPHERAL3D
template
void
computeVoronoiVolume<Dim<3>>(const FieldList<Dim<3>, Dim<3>::Vector>& position,
                             const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
                             const ConnectivityMap<Dim<3> >& connectivityMap,
                             const FieldList<Dim<3>, Dim<3>::SymTensor>& damage,
                             const std::vector<Dim<3>::FacetedVolume>& facetedBoundaries,
                             const std::vector<std::vector<Dim<3>::FacetedVolume> >& holes,
                             const std::vector<Boundary<Dim<3>>*>& boundaries,
                             const FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                             FieldList<Dim<3>, int>& surfacePoint,
                             FieldList<Dim<3>, Dim<3>::Scalar>& vol,
                             FieldList<Dim<3>, Dim<3>::Vector>& deltaMedian,
                             FieldList<Dim<3>, vector<Dim<3>::Vector>>& etaVoidPoints,
                             FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells,
                             FieldList<Dim<3>, std::vector<CellFaceFlag>>& cellFaceFlags);
#endif

}
