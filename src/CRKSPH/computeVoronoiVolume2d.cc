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
// Functor to return the mass and its derivative for use with our
// Newton-Raphson iteration.
// Note we assume here we're working in a unit radius circle (centered on the
// origin) and the gradient is aligned with and increasing in the x-direction.
//------------------------------------------------------------------------------
struct CircleMassAndGradient {
  double rho0, b;
  CircleMassAndGradient(const double rho0_in,
                        const double b_in): rho0(rho0_in),
                                            b(b_in) {}
  std::pair<double, double> operator()(const double x) const {
    CHECK(std::abs(x) <= 1.0);
    return std::make_pair(0.5*M_PI*rho0 + sqrt(1.0 - x*x)/3.0*(3.0*rho0*x + 2.0*b*(x*x - 1.0)) - rho0*acos(x),
                          2.0*(rho0 + b*x)*sqrt(1.0 - x*x));
  }
};

//------------------------------------------------------------------------------
// Integrate a linear function in a polygon.
// We do this be evaluating it for each triangle,
// using the handy relation that if f(x,y) is a linear function then the integral
// \int f(x,y) dx dy in a trianglur region is A*f(xc,yc), where A is the area of the
// triangle and (xc,yc) the triangle centroid.
// Note we implicitly use the centroid in our cell coordinates as zero.
//------------------------------------------------------------------------------
double cellIntegral(const r2d_poly& cell,
                    const double a,
                    const Dim<2>::Vector& b) {
  double result = 0.0;
  Dim<2>::Vector cent;
  int lastvert = -1, nextvert, ivert = 0, k = 0;
  while (k < cell.nverts) {
    nextvert = (cell.verts[ivert].pnbrs[0] == lastvert ?
                cell.verts[ivert].pnbrs[1] :
                cell.verts[ivert].pnbrs[0]);
    cent.x((cell.verts[ivert].pos.x + cell.verts[nextvert].pos.x)/3.0);
    cent.y((cell.verts[ivert].pos.y + cell.verts[nextvert].pos.y)/3.0);
    result += 0.5*abs(cell.verts[ivert].pos.x * cell.verts[nextvert].pos.y -
                      cell.verts[ivert].pos.y * cell.verts[nextvert].pos.x)*(a + b.dot(cent));
    lastvert = ivert;
    ivert = nextvert;
    ++k;
  }
  return result;
}

//------------------------------------------------------------------------------
// Functor to find the mass in a clipped polygon.
//------------------------------------------------------------------------------
struct PolygonClippedMassRoot {
  r2d_poly* cell0;
  mutable r2d_poly cell1;
  double rho0, targetFrac, xmin, xmax, M0;
  const Dim<2>::Vector& gradrho, searchDirection;
  mutable r2d_real mom[3];
  mutable vector<r2d_plane> clipPlane;
  PolygonClippedMassRoot(r2d_poly& celli,
                         const double rhoi,
                         const Dim<2>::Vector& gradrhoi,
                         const Dim<2>::Vector& searchDirectioni,
                         const double targetFraction):
    cell0(&celli),
    cell1(),
    rho0(rhoi),
    targetFrac(targetFraction),
    xmin(numeric_limits<double>::max()),
    xmax(-numeric_limits<double>::max()),
    M0(0.0),
    gradrho(gradrhoi),
    searchDirection(searchDirectioni),
    clipPlane(1) {

    // Find the min/max extent of the polygon along the search direction.
    double xi;
    for (unsigned j = 0; j != cell0->nverts; ++j) {
      xi = cell0->verts[j].pos.x * searchDirectioni.x() + cell0->verts[j].pos.y * searchDirection.y();
      xmin = min(xmin, xi);
      xmax = max(xmax, xi);
    }

    // Compute the total mass of the input polygon.
    M0 = cellIntegral(*cell0, rho0, gradrho);
    // mom[0] = rhoi;
    // mom[1] = gradrhoi.x();
    // mom[2] = gradrhoi.y();
    // r2d_reduce(cell0, mom, 1);
    // M0 = mom[0];
  }

  double operator()(const double x) const {
    // Clip the polygon.
    REQUIRE(x >= 0.0 and x <= 1.0);
    const double d = xmin + x*(xmax - xmin);
    if (d < 0.0) {
      clipPlane[0].n.x = -searchDirection.x();
      clipPlane[0].n.y = -searchDirection.y();
      clipPlane[0].d = d;
    } else {
      clipPlane[0].n.x = -searchDirection.x();
      clipPlane[0].n.y = -searchDirection.y();
      clipPlane[0].d = d;
    }
    cell1 = *cell0;
    r2d_clip(&cell1, &clipPlane[0], 1);
    // cout << "--------------------------------------------------------------------------------" << endl
    //      << x << endl
    //      << "Cell 0: " << endl;
    // r2d_print(cell0);
    // cout << "Clipped cell:" << endl;
    // r2d_print(&cell1);

    // Integrate the mass in the clipped cell.
    const double M1 = cellIntegral(cell1, rho0, gradrho);
    mom[0] = rho0;
    mom[1] = gradrho.x();
    mom[2] = gradrho.y();
    r2d_reduce(&cell1, mom, 1);

    // Return the difference in the mass ratio vs. the target fraction.
    // cout << " --> " << x << " " << gradrho << " " << M1 << " " << M0 << endl;
    return M1/M0 - targetFrac;
  }
};
  
//------------------------------------------------------------------------------
// Define the square distance between two r2d_vertices.
//------------------------------------------------------------------------------
double distance2(const r2d_vertex& a, const r2d_vertex& b) {
  return (FastMath::square(a.pos.x - b.pos.x) +
          FastMath::square(a.pos.y - b.pos.y));
}

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

  if (numGensGlobal > 0) {

    const Scalar rin = 0.5*kernelExtent;

    // Build an approximation of the starting kernel shape.
    const unsigned nverts = 18;
    const double dtheta = 2.0*M_PI/nverts;
    vector<Vector> verts(nverts);
    for (unsigned j = 0; j != nverts; ++j) {
      const double theta = j*dtheta;
      verts[j].x(kernelExtent*cos(theta));
      verts[j].y(kernelExtent*sin(theta));
    }
    // r2d_rvec2 verts[nverts];
    // for (unsigned j = 0; j != nverts; ++j) {
    //   const double theta = j*dtheta;
    //   verts[j].x = kernelExtent*cos(theta);
    //   verts[j].y = kernelExtent*sin(theta);
    // }
    // r2d_poly initialCell;
    // r2d_init_poly(&initialCell, verts, nverts);
    // CHECK(r2d_is_good(&initialCell));

    // Walk the points.
    r2d_real voli[1], firstmom[3];
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
        std::sort(pairPlanes.begin(), pairPlanes.end(), compareR2Dplanes);

        // Choose our seed cell shape.
        r2d_poly celli;
        if (haveBoundaries) {

          // If we have a boundary, use that for the initial cell shape.
          const vector<Facet>& facets = boundaries[nodeListi].facets();
          const unsigned nfacets = facets.size();
          r2d_rvec2 verts_bound[nfacets];
          for (unsigned j = 0; j != nfacets; ++j) {
            const Vector& vi = facets[j].point1() - ri;
            verts_bound[j].x = vi.x();
            verts_bound[j].y = vi.y();
          }
          r2d_init_poly(&celli, verts_bound, nfacets);

        } else {

          // Otherwise we use our roughly circular type.
          r2d_rvec2 verts_bound[nverts];
          for (unsigned j = 0; j != nverts; ++j) {
            const Vector vi = Hinv*verts[j];
            verts_bound[j].x = vi.x();
            verts_bound[j].y = vi.y();
          }
          r2d_init_poly(&celli, verts_bound, nverts);

        }
        CHECK2(r2d_is_good(&celli), "Bad polygon!");

        // Make a copy of the unclipped cell geometry for use in the centroidal/median filtering algorithm.
        r2d_poly doublecelli = celli;

        // Clip the local cell.
        r2d_clip(&celli, &pairPlanes[0], pairPlanes.size());

        // Check if the final polygon is entirely within our "interior" check radius.
        bool interior = true;
        {
          unsigned k = 0;
          while (interior and k != celli.nverts) {
            interior = (Hi*Vector(celli.verts[k].pos.x, celli.verts[k].pos.y)).magnitude() < rin;
            ++k;
          }
        }

        if (interior) {
          if (returnSurface) surfacePoint(nodeListi, i) = 0;

          // This is an interior point -- extract the area.
          voli[0] = 1.0;
          r2d_reduce(&celli, voli, 0);
          vol(nodeListi, i) = voli[0];

          // Apply the gradient limiter;
          gradRhoi *= phi;

          // Clip the original cell geometry again, this time using planes coincident with the neighbor points.
          for (unsigned j = 0; j != pairPlanes.size(); ++j) pairPlanes[j].d *= 2.0;
          r2d_clip(&doublecelli, &pairPlanes[0], pairPlanes.size());

          // Compute the centroidal motion.
          firstmom[0] = rhoi;
          firstmom[1] = gradRhoi.x();
          firstmom[2] = gradRhoi.y();
          r2d_reduce(&doublecelli, firstmom, 1);
          const Scalar m0 = cellIntegral(doublecelli, rhoi, gradRhoi);
          const Vector deltaCentroidi = Vector(firstmom[1], firstmom[2])/m0;


          // Is there a significant density gradient?
          if (sqrt(gradRhoi.magnitude2()*voli[0]) >= 0.025*rhoi) {

            const Vector nhat1 = gradRhoi.unitVector();
            const Vector nhat2 = Vector(-nhat1.y(), nhat1.x());
            PolygonClippedMassRoot F1(doublecelli, rhoi, gradRhoi, nhat1, 0.5);
            const double dx = F1.xmax - F1.xmin;
            const double dx1 = -F1.xmin;
            const Scalar b = gradRhoi.magnitude();
            const Scalar rho0 = rhoi - b*dx1;
            deltaMedian(nodeListi, i) = ((sqrt(2.0*rho0*rho0 + b*b*dx*dx + 2.0*b*rho0*dx)/sqrt(2.0) - rho0)*safeInvVar(b) - dx1)*nhat1 + deltaCentroidi.dot(nhat2)*nhat2;

            // // If so, we search for the median mass position within the cell.
            // // We search for the median coordinates with reference to the density gradient direction.
            // const Vector nhat1 = gradRhoi.unitVector();
            // const Vector nhat2 = Vector(-nhat1.y(), nhat1.x());
            // PolygonClippedMassRoot F1(celli, rhoi, gradRhoi, nhat1, 0.5);
            // PolygonClippedMassRoot F2(celli, rhoi, gradRhoi, nhat2, 0.5);
            // const double x1 = F1.xmin + (F1.xmax - F1.xmin)*bisectRoot(F1, 0.0, 1.0, 1.0e-5, 1.0e-5);
            // deltaMedian(nodeListi, i) = x1*nhat1 + deltaCentroidi.dot(nhat2)*nhat2;
            // // const double x2 = F2.xmin + (F2.xmax - F2.xmin)*bisectRoot(F2, 0.0, 1.0, 1.0e-5, 1.0e-5);
            // // deltaMedian(nodeListi, i) = x1*nhat1 + x2*nhat2;
            // // cout << " **> " << i << " " << deltaMedian(nodeListi, i) << " " << phi << " " << gradRhoi << " " << x1 << " " << x2 << " " << (x1 - F1.xmin)/(F1.xmax - F1.xmin) << " " << (x2 - F2.xmin)/(F2.xmax - F2.xmin) << endl;
            // // cout << "================================================================================" << endl;

          } else {

            // Otherwise just use the centroid.
            deltaMedian(nodeListi, i) = deltaCentroidi;
            // cout << " CENTROIDAL**> " << i << " " << deltaMedian(nodeListi, i) << " " << phi << " " << gradRhoi << endl;

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

        // If requested, we can return the cell geometries.
        if (returnCells) {
          vector<Vector> verts;
          verts.reserve(celli.nverts);
          vector<vector<unsigned> > facetIndices; // (celli.nverts, vector<unsigned>(2));
          int lastvert = -1, nextvert, ivert = 0, j = 0, k = 0;
          const Scalar tol = 1.0e-8*sqrt(Hdeti);
          while (k < celli.nverts) {
            if (lastvert == -1 or
                distance2(celli.verts[ivert], celli.verts[lastvert]) > tol) {
              verts.push_back(Vector(celli.verts[ivert].pos.x, celli.verts[ivert].pos.y) + ri);
              facetIndices.push_back(vector<unsigned>(2));
              CHECK(facetIndices.size() == j + 1);
              facetIndices[j][0] = j;
              facetIndices[j][1] = j + 1;
              ++j;
            }
            nextvert = (celli.verts[ivert].pnbrs[0] == lastvert ?
                        celli.verts[ivert].pnbrs[1] :
                        celli.verts[ivert].pnbrs[0]);
            lastvert = ivert;
            ivert = nextvert;
            ++k;
          }
          if ((verts.back() - verts.front()).magnitude2() <= tol) { // The first and last may have snuck in degenerate.
            verts.pop_back();
            facetIndices.pop_back();
          }
          facetIndices.back()[1] = 0;
          // std::copy(verts.begin(), verts.end(), std::ostream_iterator<Dim<2>::Vector>(std::cout, " "));
          // std::cout << endl;
          cells(nodeListi, i) = FacetedVolume(verts, facetIndices);
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
