//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
extern "C" {
#include "r3d/r2d.h"
}

#include <algorithm>
#include <utility>

#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/FastMath.hh"
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
// A special comparator to sort r2d planes by distance.
//------------------------------------------------------------------------------
bool compareR2Dplanes(const r2d_plane& lhs, const r2d_plane& rhs) {
  return lhs.d < rhs.d;
}

//------------------------------------------------------------------------------
// Find the 1D extent of an R2D cell along the given direction.
//------------------------------------------------------------------------------
void findPolygonExtent(double& xmin, double& xmax, const Dim<2>::Vector& nhat, const r2d_poly& celli) {
  REQUIRE(fuzzyEqual(nhat.magnitude(), 1.0));
  double xi;
  xmin = 0.0;
  xmax = 0.0;
  for (unsigned i = 0; i != celli.nverts; ++i) {
    xi = (celli.verts[i].pos.x * nhat.x() +
          celli.verts[i].pos.y * nhat.y());
    xmin = std::min(xmin, xi);
    xmax = std::max(xmax, xi);
  }
  xmin = std::min(0.0, xmin);
  xmax = std::max(0.0, xmax);
}

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

}           // anonymous namespace

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
void
computeVoronoiVolume(const FieldList<Dim<2>, Dim<2>::Vector>& position,
                     const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
                     const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>& rho,
                     const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& gradRho,
                     const ConnectivityMap<Dim<2> >& connectivityMap,
                     const Dim<2>::Scalar kernelExtent,
                     const std::vector<Dim<2>::FacetedVolume>& boundaries,
                     const std::vector<std::vector<Dim<2>::FacetedVolume> >& holes,
                     FieldList<Dim<2>, int>& surfacePoint,
                     FieldList<Dim<2>, Dim<2>::Scalar>& vol,
                     FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& deltaMedian,
                     FieldSpace::FieldList<Dim<2>, Dim<2>::FacetedVolume>& cells) {

  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;
  typedef Dim<2>::FacetedVolume::Facet Facet;

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();
  const unsigned numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());
  const unsigned numBounds = boundaries.size();
  const bool haveBoundaries = numBounds == numNodeLists;
  const bool returnSurface = surfacePoint.size() == numNodeLists;
  const bool returnCells = cells.size() == numNodeLists;

  REQUIRE(numBounds == 0 or numBounds == numNodeLists);
  REQUIRE(holes.size() == numBounds);

  if (numGensGlobal > 0) {

    const Scalar rin2 = 0.25*kernelExtent*kernelExtent;

    // Build an approximation of the starting kernel shape.
    const unsigned nverts = 18;
    const double dtheta = 2.0*M_PI/nverts;
    r2d_rvec2 verts[nverts];
    for (unsigned j = 0; j != nverts; ++j) {
      const double theta = j*dtheta;
      verts[j].x = kernelExtent*cos(theta);
      verts[j].y = kernelExtent*sin(theta);
    }
    r2d_poly initialCell;
    r2d_init_poly(&initialCell, verts, nverts);
    CHECK(r2d_is_good(&initialCell));

    // Walk the points.
    r2d_real firstmom[3];
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = vol[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {

        // const bool barf = (i == 1);

        const Vector& ri = position(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const Scalar rhoi = rho(nodeListi, i);
        const Vector& gradRhoi = gradRho(nodeListi, i);
        const Vector grhat = gradRhoi.unitVector();
        const Scalar Hdeti = Hi.Determinant();
        const SymTensor Hinv = Hi.Inverse();

        // if (barf) cerr << " --> " << i << " " << ri << endl;

        // Grab this points neighbors and build all the planes.
        // We simultaneously build a very conservative limiter for the density gradient.
        Scalar phi = 1.0;
        vector<r2d_plane> pairPlanes;
        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          for (vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
               jItr != fullConnectivity[nodeListj].end();
               ++jItr) {
            const unsigned j = *jItr;
            const Vector& rj = position(nodeListj, j);
            const Scalar rhoj = rho(nodeListj, j);

            // Build the planes for our clipping half-spaces.
            const Vector rij = ri - rj;
            const Vector nhat = rij.unitVector();
            pairPlanes.push_back(r2d_plane());
            pairPlanes.back().n.x = nhat.x();
            pairPlanes.back().n.y = nhat.y();
            pairPlanes.back().d = 0.5*rij.magnitude();

            // Check the density gradient limiter.
            const Scalar fdir = FastMath::pow4(rij.unitVector().dot(grhat));
            // phi = min(phi, max(0.0, max(1.0 - fdir, rij.dot(gradRhoi)*safeInv(rhoi - rhoj))));
          }
        }

        // If provided boundaries, we implement them as additional neighbor clipping planes.
        if (haveBoundaries) {
          const vector<Facet>& facets = boundaries[nodeListi].facets();
          for (const Facet& facet: facets) {
            const Vector p = facet.closestPoint(ri);
            Vector rij = ri - p;
            if (rij.magnitude2() < kernelExtent*kernelExtent) {
              Vector nhat;
              if (rij.magnitude() < 1.0e-3*facet.area()) {
                rij.Zero();
                nhat = -facet.normal();
              } else {
                nhat = rij.unitVector();
              }
              pairPlanes.push_back(r2d_plane());
              pairPlanes.back().n.x = nhat.x();
              pairPlanes.back().n.y = nhat.y();
              pairPlanes.back().d = rij.magnitude();
            }
          }

          // Same thing with holes.
          for (const FacetedVolume& hole: holes[nodeListi]) {
            const vector<Facet>& facets = hole.facets();
            for (const Facet& facet: facets) {
              const Vector p = facet.closestPoint(ri);
              Vector rij = ri - p;
              if (rij.magnitude2() < kernelExtent*kernelExtent) {
                Vector nhat;
                if (rij.magnitude2() < 1.0e-3*facet.area()) {
                  rij.Zero();
                  nhat = facet.normal();
                } else {
                  nhat = rij.unitVector();
                }
                pairPlanes.push_back(r2d_plane());
                pairPlanes.back().n.x = nhat.x();
                pairPlanes.back().n.y = nhat.y();
                pairPlanes.back().d = rij.magnitude();
              }
            }
          }
        }

        // Sort the planes by distance -- let's us clip more efficiently.
        std::sort(pairPlanes.begin(), pairPlanes.end(), compareR2Dplanes);

        // Initialize our seed cell shape.
        r2d_poly celli = initialCell;
        for (r2d_vertex& vert: celli.verts) {
          const Vector vi = Hinv*Vector(vert.pos.x, vert.pos.y);
          vert.pos.x = vi.x();
          vert.pos.y = vi.y();
        }
        CHECK2(r2d_is_good(&celli), "Bad initial polygon!");

        // Clip the local cell.
        r2d_clip(&celli, &pairPlanes[0], pairPlanes.size());
        CHECK(celli.nverts > 0);

        // if (barf) r2d_print(&celli);

        // Check if the final polygon is entirely within our "interior" check radius.
        bool interior = true;
        {
          unsigned k = 0;
          while (interior and k != celli.nverts) {
            interior = (Hi*Vector(celli.verts[k].pos.x, celli.verts[k].pos.y)).magnitude2() < rin2;
            ++k;
          }
        }

        if (interior) {
          if (returnSurface) surfacePoint(nodeListi, i) = 0;

          // Compute the centroidal motion and area.
          r2d_reduce(&celli, firstmom, 1);
          CHECK(firstmom[0] > 0.0);
          vol(nodeListi, i) = firstmom[0];
          const Vector deltaCentroidi = Vector(firstmom[1], firstmom[2])/firstmom[0];
          // if (barf) cerr << "     " << deltaCentroidi << " " << ri + deltaCentroidi << endl;

          // // Apply the gradient limiter;
          // gradRhoi *= phi;

          // Is there a significant density gradient?
          if (sqrt(gradRhoi.magnitude2()*vol(nodeListi, i)) >= 1e-8*rhoi) {

            const Vector nhat1 = gradRhoi.unitVector();
            double x1, x2;
            findPolygonExtent(x1, x2, nhat1, celli);
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

            // // This version simply tries rho^2 weighting.
            // deltaMedian(nodeListi, i) = ((0.5*rhoi*(x2*x2 - x1*x1) +
            //                               2.0/3.0*rhoi*b*(x2*x2*x2 - x1*x1*x1) +
            //                               0.25*b*b*(x2*x2*x2*x2 - x1*x1*x1*x1))/
            //                              (pow3(rhoi + b*x2) - pow3(rhoi + b*x1)/(3.0*b)))*nhat1 - deltaCentroidi.dot(nhat1)*nhat1 + deltaCentroidi;

          } else {

            // Otherwise just use the centroid.
            deltaMedian(nodeListi, i) = deltaCentroidi;

          }

          // OK, this is an interior point from the perspective that it was clipped within our critical
          // radius on all sides.  However, if we have a bounding polygon we may still want to call it a
          // surface if in fact there are still facets from that bounding polygon on this cell.
          if (haveBoundaries and returnSurface) {
            unsigned j = 0;
            while (interior and j != celli.nverts) {
              interior = not pointOnPolygon(ri + Vector(celli.verts[j].pos.x, celli.verts[j].pos.y),
                                            boundaries[nodeListi].vertices(),
                                            1.0e-8);
              ++j;
            }

            if (not interior) {
              // This is a point that touches the bounding polygon.  Flag it as surface.
              surfacePoint(nodeListi, i) = 1;
            }
          }

        } else {

          // This point touches a free boundary, so flag it.
          if (returnSurface) surfacePoint(nodeListi, i) = 1;
          deltaMedian(nodeListi, i) = Vector::zero;

        }

        // Check if the candidate motion is still in the boundary.  If not, project back.
        if (haveBoundaries) {
          if (not boundaries[nodeListi].contains(ri + deltaMedian(nodeListi, i))) {
            deltaMedian(nodeListi, i) = boundaries[nodeListi].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
          }
          for (unsigned ihole = 0; ihole != holes[nodeListi].size(); ++ihole) {
            if (holes[nodeListi][ihole].contains(ri + deltaMedian(nodeListi, i))) {
              deltaMedian(nodeListi, i) = holes[nodeListi][ihole].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
            }
          }
        }

        // If requested, we can return the cell geometries.
        if (returnCells) {
          r2d_poly_to_polygon(celli, 1.0e-8/sqrt(Hdeti), cells(nodeListi, i));
          cells(nodeListi, i) += ri;
        }
      }
    }

    // // Lastly, for any points labeled suface we concoct a sampled median motion by simply interpolating from the surrounding points.
    // // We'll use a really simple low-order Shepard's function for this.
    // if (numBounds == 0) {
    //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    //     const unsigned n = vol[nodeListi]->numInternalElements();
    //     const Neighbor<Dim<2> >& neighbor = position[nodeListi]->nodeListPtr()->neighbor();
    //     for (unsigned i = 0; i != n; ++i) {
    //       const Vector& ri = position(nodeListi, i);
    //       Scalar wsumi = 0.0;
    //       const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
    //       for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
    //         for (vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
    //              jItr != fullConnectivity[nodeListj].end();
    //              ++jItr) {
    //           const unsigned j = *jItr;
    //           if (surfacePoint(nodeListj, j) == 0) {
    //             const Vector& rj = position(nodeListj, j);
    //             const Scalar wi = safeInvVar((ri - rj).magnitude2());
    //             wsumi += wi;
    //             deltaMedian(nodeListi, i) += wi*deltaMedian(nodeListj, j);
    //           }
    //         }
    //       }
    //       deltaMedian(nodeListi, i) *= safeInvVar(wsumi);
    //     }
    //   }
    // }

  }
}

}
}
