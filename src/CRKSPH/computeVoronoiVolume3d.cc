//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
extern "C" {
#include "r3d/r3d.h"
}

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
#include "Utilities/r3d_utils.hh"

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
// A special comparator to sort r3d planes by distance.
//------------------------------------------------------------------------------
bool compareR3Dplanes(const r3d_plane& lhs, const r3d_plane& rhs) {
  return lhs.d < rhs.d;
}

//------------------------------------------------------------------------------
// Find the 1D extent of an R3D cell along the given direction.
//------------------------------------------------------------------------------
void findPolyhedronExtent(double& xmin, double& xmax, const Dim<3>::Vector& nhat, const r3d_poly& celli) {
  REQUIRE(fuzzyEqual(nhat.magnitude(), 1.0));
  double xi;
  xmin = 0.0;
  xmax = 0.0;
  for (unsigned i = 0; i != celli.nverts; ++i) {
    xi = (celli.verts[i].pos.x * nhat.x() +
          celli.verts[i].pos.y * nhat.y() +
          celli.verts[i].pos.z * nhat.z());
    xmin = std::min(xmin, xi);
    xmax = std::max(xmax, xi);
  }
  xmin = std::min(0.0, xmin);
  xmax = std::max(0.0, xmax);
}

// //------------------------------------------------------------------------------
// // Return a Spheral GeomPolyhedron from an R3D polyhedron.
// //------------------------------------------------------------------------------
// Dim<3>::FacetedVolume
// r3d_poly_to_polyhedron(const r3d_poly& celli,
//                        const Dim<3>::Vector& offset,
//                        const double tol) {

//   using std::vector;
//   typedef Dim<3>::Scalar Scalar;
//   typedef Dim<3>::Vector Vector;
//   typedef Dim<3>::SymTensor SymTensor;
//   typedef Dim<3>::FacetedVolume FacetedVolume;
//   typedef Dim<3>::FacetedVolume::Facet Facet;

//   // We're going to cheat here a bit.  We know that currently computeVoronoiVolume3d only returns
//   // convex polyhedra, so we'll take the easy way out and just build the convex hull of the
//   // vertices.
//   vector<Vector> verts;
//   for (auto i = 0; i != celli.nverts; ++i) verts.push_back(Vector(celli.verts[i].pos.x,
//                                                                   celli.verts[i].pos.y,
//                                                                   celli.verts[i].pos.z));
//   return FacetedVolume(verts);
// }

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
                     const Dim<3>::Scalar kernelExtent,
                     const std::vector<Dim<3>::FacetedVolume>& boundaries,
                     const std::vector<std::vector<Dim<3>::FacetedVolume> >& holes,
                     const FieldSpace::FieldList<Dim<3>, Dim<3>::Scalar>& weight,
                     FieldList<Dim<3>, int>& surfacePoint,
                     FieldList<Dim<3>, Dim<3>::Scalar>& vol,
                     FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& deltaMedian,
                     FieldSpace::FieldList<Dim<3>, Dim<3>::FacetedVolume>& cells) {

  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::FacetedVolume FacetedVolume;
  typedef Dim<3>::FacetedVolume::Facet Facet;

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();
  const unsigned numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());
  const unsigned numBounds = boundaries.size();
  const bool haveBoundaries = numBounds == numNodeLists;
  const bool haveWeights = weight.size() == numNodeLists;
  const bool returnSurface = surfacePoint.size() == numNodeLists;
  const bool returnCells = cells.size() == numNodeLists;

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

  if (numGensGlobal > 0) {

    // (Square) of the distance to a facet in an icosahedon.
    const Scalar rin2 = 0.25*kernelExtent*kernelExtent * (15.0*(5.0+sqrt(5.0)))/900.0;

    // Build an approximation of the starting kernel shape (in eta space) as an icosahedron.
    const unsigned nverts = 12;
    const unsigned nfaces = 20;
    r3d_int faces[nfaces][3] = {
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
      {9, 8, 1},
    };
    r3d_int** facesp = new r3d_int*[nfaces];
    for (unsigned j = 0; j != nfaces; ++j) {
      facesp[j] = new r3d_int[3];
      for (unsigned k = 0; k != 3; ++k) facesp[j][k] = faces[j][k];
    }
    r3d_int nvertsperface[nfaces] = {  // Array of number of vertices per face.
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
    };
    const double t = (1.0 + sqrt(5.0)) / 2.0;
    r3d_rvec3 verts[nverts];           // Array of vertex coordinates.
    verts[0].x =  -1; verts[0].y =  t; verts[0].z =   0;
    verts[1].x =   1; verts[1].y =  t; verts[1].z =   0;
    verts[2].x =  -1; verts[2].y = -t; verts[2].z =   0;
    verts[3].x =   1; verts[3].y = -t; verts[3].z =   0;
    verts[4].x =   0; verts[4].y = -1; verts[4].z =   t;
    verts[5].x =   0; verts[5].y =  1; verts[5].z =   t;
    verts[6].x =   0; verts[6].y = -1; verts[6].z =  -t;
    verts[7].x =   0; verts[7].y =  1; verts[7].z =  -t;
    verts[8].x =   t; verts[8].y =  0; verts[8].z =  -1;
    verts[9].x =   t; verts[9].y =  0; verts[9].z =   1;
    verts[10].x = -t; verts[10].y = 0; verts[10].z = -1;
    verts[11].x = -t; verts[11].y = 0; verts[11].z =  1;
    r3d_poly initialCell;
    r3d_init_poly(&initialCell, verts, nverts, facesp, nvertsperface, nfaces);
    CHECK(r3d_is_good(&initialCell));

    // Deallocate that damn memory. I hate this syntax, but don't know enough C to know if there's a better way.
    for (unsigned j = 0; j != nfaces; ++j) delete[] facesp[j];
    delete[] facesp;

    // Scale the template icosahedron to have the initial volume of a sphere of radius kernelExtent.
    r3d_real voli[1], firstmom[4];
    r3d_reduce(&initialCell, voli, 0);
    CHECK(voli[0] > 0.0);
    const double volscale = Dim<3>::rootnu(4.0/3.0*M_PI/voli[0])*kernelExtent;
    r3d_scale(&initialCell, volscale);
    CHECK(r3d_is_good(&initialCell));
    BEGIN_CONTRACT_SCOPE
    {
      r3d_reduce(&initialCell, voli, 0);
      CHECK2(fuzzyEqual(voli[0], 4.0/3.0*M_PI*Dim<3>::pownu(kernelExtent), 1.0e-10), voli[0] << " " << 4.0/3.0*M_PI*Dim<3>::pownu(kernelExtent) << " " << volscale);
    }
    END_CONTRACT_SCOPE
    
    // Walk the points.
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = vol[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {

        const Vector& ri = position(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const Scalar rhoi = rho(nodeListi, i);
        Vector gradRhoi = gradRho(nodeListi, i);
        const Vector grhat = gradRhoi.unitVector();
        const Scalar Hdeti = Hi.Determinant();
        const SymTensor Hinv = Hi.Inverse();
        const Scalar weighti = haveWeights ? weight(nodeListi, i) : 1.0;

        // t0 = std::clock();

        // Grab this points neighbors and build all the planes.
        // We simultaneously build a very conservative limiter for the density gradient.
        // Scalar phi = 1.0;
        vector<r3d_plane> pairPlanes;
        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          for (vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
               jItr != fullConnectivity[nodeListj].end();
               ++jItr) {
            const unsigned j = *jItr;
            const Vector& rj = position(nodeListj, j);
            const Scalar rhoj = rho(nodeListj, j);
            const Scalar weightj = haveWeights ? weight(nodeListj, j) : 1.0;

            // Build the planes for our clipping half-spaces.
            const Vector rij = ri - rj;
            const Vector nhat = rij.unitVector();
            const Scalar wij = weighti/(weightj + weighti);
            pairPlanes.push_back(r3d_plane());
            pairPlanes.back().n.x = nhat.x();
            pairPlanes.back().n.y = nhat.y();
            pairPlanes.back().n.z = nhat.z();
            pairPlanes.back().d = wij*rij.magnitude();

            // // Check the density gradient limiter.
            // const Scalar fdir = FastMath::pow4(rij.unitVector().dot(grhat));
            // phi = min(phi, max(0.0, max(1.0 - fdir, rij.dot(gradRhoi)*safeInv(rhoi - rhoj))));
          }
        }

        // tplanesneighbors += std::clock() - t0;
        // t0 = std::clock();

        // If provided holes, we implement them as additional neighbor clipping planes.
        if (haveBoundaries) {
          const vector<Facet>& facets = boundaries[nodeListi].facets();
          CHECK(boundaries[nodeListi].contains(ri, false));
          for (const Facet& facet: facets) {
            const Vector p = facet.closestPoint(ri);
            Vector rij = ri - p;
            if ((Hi*rij).magnitude2() < kernelExtent*kernelExtent) {
              Vector nhat;
              if (rij.magnitude() < 1.0e-5*facet.area()) {
                rij.Zero();
                nhat = -facet.normal();
              } else {
                nhat = rij.unitVector();
              }
              pairPlanes.push_back(r3d_plane());
              pairPlanes.back().n.x = nhat.x();
              pairPlanes.back().n.y = nhat.y();
              pairPlanes.back().n.z = nhat.z();
              pairPlanes.back().d = rij.magnitude();
            }
          }

          // Same thing with holes.
          for (const FacetedVolume& hole: holes[nodeListi]) {
            CHECK(not hole.contains(ri));
            const vector<Facet>& facets = hole.facets();
            for (const Facet& facet: facets) {
              const Vector p = facet.closestPoint(ri);
              Vector rij = ri - p;
              if ((Hi*rij).magnitude2() < kernelExtent*kernelExtent) {
                Vector nhat;
                if (rij.magnitude2() < 1.0e-5*facet.area()) {
                  rij.Zero();
                  nhat = facet.normal();
                } else {
                  nhat = rij.unitVector();
                }
                pairPlanes.push_back(r3d_plane());
                pairPlanes.back().n.x = nhat.x();
                pairPlanes.back().n.y = nhat.y();
                pairPlanes.back().n.z = nhat.z();
                pairPlanes.back().d = rij.magnitude();
              }
            }
          }
        }

        // tplanesboundaries += std::clock() - t0;
        // t0 = std::clock();

        // Sort the planes by distance -- let's us clip more efficiently.
        std::sort(pairPlanes.begin(), pairPlanes.end(), compareR3Dplanes);

        // tplanesort += std::clock() - t0;
        // t0 = std::clock();

        // Initialize our seed cell shape.  If we have a boundary use that, otherwise the nominal kernel boundary.
        r3d_poly celli;
        if (false) { // (haveBoundaries) {
          polyhedron_to_r3d_poly(boundaries[nodeListi] - ri, celli);
        } else {
          celli = initialCell;
          for (unsigned k = 0; k != celli.nverts; ++k) {
            r3d_vertex& vert = celli.verts[k];
            const Vector vi = Hinv*Vector(vert.pos.x, vert.pos.y, vert.pos.z);
            vert.pos.x = vi.x();
            vert.pos.y = vi.y();
            vert.pos.z = vi.z();
          }
        }
        CHECK2(r3d_is_good(&celli), "Bad initial polyhedron!");

        // Clip the local cell.
        r3d_clip(&celli, &pairPlanes[0], pairPlanes.size());
        CHECK2(r3d_is_good(&celli), "Bad polyhedron!");
        r3d_reduce(&celli, firstmom, 1);
        CHECK(firstmom[0] > 0.0);
        vol(nodeListi, i) = firstmom[0];
        const Vector deltaCentroidi = Vector(firstmom[1], firstmom[2], firstmom[3])/firstmom[0];
        deltaMedian(nodeListi, i) = deltaCentroidi;

        // tclip += std::clock() - t0;

        // Check if the final polyhedron is entirely within our "interior" check radius.
        bool interior = true;
        // t0 = std::clock();
        {
          unsigned k = 0;
          while (interior and k != celli.nverts) {
            interior = (Hi*Vector(celli.verts[k].pos.x, celli.verts[k].pos.y, celli.verts[k].pos.z)).magnitude2() < rin2;
            ++k;
          }
        }
        // tinterior += std::clock() - t0;

        if (interior) {
          // t0 = std::clock();
          if (returnSurface) surfacePoint(nodeListi, i) = 0;

          // // Apply the gradient limiter;
          // gradRhoi *= phi;

          // Is there a significant density gradient?
          if (gradRhoi.magnitude()*Dim<3>::rootnu(vol(nodeListi, i)) >= 1e-8*rhoi) {

            const Vector nhat1 = gradRhoi.unitVector();
            double x1, x2;
            findPolyhedronExtent(x1, x2, nhat1, celli);
            CHECK2(x1 <= 0.0 and x2 >= 0.0, nodeListi << " " << i << " " << ri << " " << x1 << " " << x2);
            const Scalar b = gradRhoi.magnitude();

            // This version uses the medial position.
            const Scalar thpt = sqrt(abs(rhoi*rhoi + rhoi*b*(x1 + x2) + 0.5*b*b*(x1*x1 + x2*x2)));
            const Scalar xm1 = -(rhoi + thpt)/b;
            const Scalar xm2 = (-rhoi + thpt)/b;
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
            unsigned j = 0;
            while (interior and j != celli.nverts) {
              interior = not pointOnPolyhedron(ri + Vector(celli.verts[j].pos.x, celli.verts[j].pos.y, celli.verts[j].pos.z),
                                               boundaries[nodeListi],
                                               1.0e-8);
              ++j;
            }

            if (not interior) {
              // This is a point that touches the bounding polygon.  Flag it as surface.
              if (returnSurface) surfacePoint(nodeListi, i) = 1;
            }
          }
          // tsurface += std::clock() - t0;

        } else {

          // This point touches a free boundary, so flag it.
          if (returnSurface) surfacePoint(nodeListi, i) = 1;

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

        // If requested, we can return the cell geometries.
        if (returnCells) {
          // t0 = std::clock();
          r3d_poly_to_polyhedron(celli, 1.0e-20/max(1.0, Dim<3>::rootnu(Hdeti)), cells(nodeListi, i));
          cells(nodeListi, i) += ri;
          // tcell += std::clock() - t0;
        }
      }
    }
  }

  // ttotal = std::clock() - ttotal;
  // if (Process::getRank() == 0) cout << "computeVoronoiVolume2d timing: " 
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
}
    
}
}
