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

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

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
  const unsigned nverts = celli.nverts;
  double xi;
  xmin = std::numeric_limits<double>::max();
  xmax = -std::numeric_limits<double>::max();
  for (unsigned i = 0; i != nverts; ++i) {
    xi = (celli.verts[i].pos.x * nhat.x() +
          celli.verts[i].pos.y * nhat.y());
    xmin = std::min(xmin, xi);
    xmax = std::max(xmax, xi);
  }
}

// //------------------------------------------------------------------------------
// // Convert a Spheral polygon to an R2D polygon.  Since our polygon may be
// // non-trivial (include holes) we can't just use r2d_init_poly.
// //------------------------------------------------------------------------------
// void Polygon_to_r2d_poly(r2d_poly& cell_r2d,
//                          const Dim<2>::FacetedVolume& cell_spheral,
//                          const vector<Dim<2>::FacetedVolume>& holes_spheral) {
//   typedef Dim<2>::Vector Vector;
//   typedef Dim<2>::FacetedVolume::Facet Facet;

//   // Build the outer boundary.
//   const vector<Vector>& vertices = cell_spheral.vertices();
//   const vector<Facet>& facets = cell_spheral.facets();
//   const unsigned nverts = vertices.size();
//   CHECK(nverts <= R2D_MAX_VERTS);
//   CHECK(facets.size() == nverts);
//   cell_r2d.nverts = nverts;
//   for (unsigned i = 0; i != nverts; ++i) {
//     cell_r2d.verts[i].pos.x = vertices[i].x();
//     cell_r2d.verts[i].pos.y = vertices[i].y();
//     cell_r2d.verts[facets[i].ipoint1()].pnbrs[0] = facets[i].ipoint2();
//     cell_r2d.verts[facets[i].ipoint2()].pnbrs[1] = facets[i].ipoint1();
//   }

//   // Add any holes by reversing their geometry order to be CW.
//   const unsigned nholes = holes_spheral.size();
//   for (unsigned ihole = 0; ihole != nholes; ++ihole) {
//     const vector<Vector>& vertices = holes_spheral[ihole].vertices();
//     const vector<Facet>& facets = holes_spheral[ihole].facets();
//     const unsigned nverts_old = cell_r2d.nverts;
//     const unsigned nverts = vertices.size();
//     CHECK(nverts_old + nverts <= R2D_MAX_VERTS);
//     CHECK(facets.size() == nverts);
//     cell_r2d.nverts = nverts_old + nverts;
//     for (unsigned i = 0; i != nverts; ++i) {
//       cell_r2d.verts[nverts_old + i].pos.x = vertices[i].x();
//       cell_r2d.verts[nverts_old + i].pos.y = vertices[i].y();
//       cell_r2d.verts[nverts_old + facets[i].ipoint1()].pnbrs[0] = nverts_old + facets[i].ipoint2();
//       cell_r2d.verts[nverts_old + facets[i].ipoint2()].pnbrs[1] = nverts_old + facets[i].ipoint1();
//     }
//   }

//   ENSURE2(r2d_is_good(&cell_r2d), "Bad polygon transformation.");
// }

// //------------------------------------------------------------------------------
// // Integrate a linear function in a polygon.
// // We do this be evaluating it for each triangle,
// // using the handy relation that if f(x,y) is a linear function then the integral
// // \int f(x,y) dx dy in a trianglur region is A*f(xc,yc), where A is the area of the
// // triangle and (xc,yc) the triangle centroid.
// // Note we implicitly use the centroid in our cell coordinates as zero.
// //------------------------------------------------------------------------------
// double cellIntegral(const r2d_poly& cell,
//                     const double a,
//                     const Dim<2>::Vector& b) {
//   double result = 0.0;
//   Dim<2>::Vector cent;
//   int lastvert = -1, nextvert, ivert = 0, k = 0;
//   while (k < cell.nverts) {
//     nextvert = (cell.verts[ivert].pnbrs[0] == lastvert ?
//                 cell.verts[ivert].pnbrs[1] :
//                 cell.verts[ivert].pnbrs[0]);
//     cent.x((cell.verts[ivert].pos.x + cell.verts[nextvert].pos.x)/3.0);
//     cent.y((cell.verts[ivert].pos.y + cell.verts[nextvert].pos.y)/3.0);
//     result += 0.5*abs(cell.verts[ivert].pos.x * cell.verts[nextvert].pos.y -
//                       cell.verts[ivert].pos.y * cell.verts[nextvert].pos.x)*(a + b.dot(cent));
//     lastvert = ivert;
//     ivert = nextvert;
//     ++k;
//   }
//   return result;
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
        Vector gradRhoi = gradRho(nodeListi, i);
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
            phi = min(phi, max(0.0, max(1.0 - fdir, rij.dot(gradRhoi)*safeInv(rhoi - rhoj))));
          }
        }

        // If provided boundaries, we implement them as additional neighbor clipping planes.
        if (haveBoundaries) {
          const vector<Facet>& facets = boundaries[nodeListi].facets();
          BOOST_FOREACH(const Facet& facet, facets) {
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
          BOOST_FOREACH(const FacetedVolume& hole, holes[nodeListi]) {
            const vector<Facet>& facets = hole.facets();
            BOOST_FOREACH(const Facet& facet, facets) {
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
        BOOST_FOREACH(r2d_vertex& vert, celli.verts) {
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

          // Apply the gradient limiter;
          gradRhoi *= phi;

          // Is there a significant density gradient?
          if (sqrt(gradRhoi.magnitude2()*vol(nodeListi, i)) >= 0.025*rhoi) {

            const Vector nhat1 = gradRhoi.unitVector();
            double dx1, dx2;
            findPolygonExtent(dx1, dx2, nhat1, celli);
            dx1 = -dx1;
            CHECK(dx1 >= 0. and dx2 >= 0.0);
            const Scalar b = gradRhoi.magnitude();
            deltaMedian(nodeListi, i) = (sqrt(abs(rhoi*rhoi + b*rhoi*(dx2 - dx1) + b*b*(dx1*dx1 + dx2*dx2))) - rhoi)/b*nhat1 -  deltaCentroidi.dot(nhat1)*nhat1 + deltaCentroidi;

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
        // Note, R2D leaves lots of degeneracies in the cell points/edges, so we do this in two passes.  First,
        // read all the vertices in CCW order and build a linked list pointing to the next one.  Then we
        // go over these points and remove any degeneracies by updating just the linked list to loop over
        // unique vertices.
        if (returnCells) {

          // if (barf) { // BLAGO
          //   cerr << "Raw verts: " << endl;
          //   for (unsigned j = 0; j != celli.nverts; ++j) {
          //     cerr << " --> " << celli.verts[j].pos.x + ri.x() << " " << celli.verts[j].pos.y + ri.y() << endl;
          //   }
          // } // BLAGO

          // Read out the R2D cell in CCW order.  We have to scan for the positive loop of edges though.
          vector<Vector> verts;
          vector<int> vertcheck(celli.nverts, 0);
          {
            int nextvert, ivert, firstvert;
            double area = -1.0;
            while (area < 0.0) {
              area = 0.0;

              // Find the first unused vertex.
              firstvert = 0;
              while (firstvert != celli.nverts and vertcheck[firstvert] == 1) firstvert++;
              CHECK(firstvert != celli.nverts);

              // Read out the loop of vertices.
              ivert = firstvert;
              nextvert = -1;
              verts.clear();
              while (nextvert != firstvert) {
                verts.push_back(Vector(celli.verts[ivert].pos.x,
                                       celli.verts[ivert].pos.y));
                vertcheck[ivert] = 1;
                nextvert = celli.verts[ivert].pnbrs[0];
                // if (barf) cerr << " **> " << (verts.back() + ri) << " " << ivert << " "  << nextvert << endl;
                area += (celli.verts[ivert].pos.x * celli.verts[nextvert].pos.y -
                         celli.verts[ivert].pos.y * celli.verts[nextvert].pos.x);
                ivert = nextvert;
              }
            }
            // if (barf) cerr << " area : " << area << endl;
          }

          // Flag any redundant vertices to not be used.
          vector<int> usevert(verts.size(), 1);
          const Scalar tol = 1.0e-8/sqrt(Hdeti);
          for (int j = 0; j != verts.size() - 1; ++j) {
            for (int k = j + 1; k != verts.size(); ++k) {
              if (usevert[k] == 1 and (verts[j] - verts[k]).magnitude2() < tol) usevert[k] = 0;
            }
          }

          // Now we can read out the vertices we're actually using and build the return polygon.
          vector<Vector> uniqueVerts;
          vector<vector<unsigned> > facetIndices;
          int k = 0;
          for (int j = 0; j != verts.size(); ++j) {
            if (usevert[j] == 1) {
              uniqueVerts.push_back(ri + verts[j]);
              facetIndices.push_back(vector<unsigned>(2));
              facetIndices.back()[0] = k;
              facetIndices.back()[1] = ++k;
            }
          }
          facetIndices.back()[1] = 0;
          CHECK(uniqueVerts.size() >= 3);

          // // Check the dang things are in CCW order.
          // double area = 0.0;
          // for (int j = 0; j != uniqueVerts.size(); ++j) {
          //   area += ((uniqueVerts[facetIndices[j][0]] - ri).cross(uniqueVerts[facetIndices[j][1]] - ri)).z();
          // }
          // if (area < 0.0) std::reverse(uniqueVerts.begin(), uniqueVerts.end());

          // if (barf) {
          //   cout << " --> " << i << " : ";
          //   std::copy(verts.begin(), verts.end(), std::ostream_iterator<Dim<2>::Vector>(std::cout, " "));
          //   std::cout << endl;
          //   cout << " --> " << i << " : ";
          //   std::copy(uniqueVerts.begin(), uniqueVerts.end(), std::ostream_iterator<Dim<2>::Vector>(std::cout, " "));
          //   std::cout << endl;
          // }
          cells(nodeListi, i) = FacetedVolume(uniqueVerts, facetIndices);
          // if (barf) cerr << cells(nodeListi, i) << endl;
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
