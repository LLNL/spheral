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
#include "Utilities/newtonRaphson.hh"

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
  const bool returnCells = cells.size() == numNodeLists;

  REQUIRE(numBounds == 0 or numBounds == numNodeLists);

  if (numGensGlobal > 0) {

    // Start out assuming all points are internal.
    surfacePoint = 0;
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
      const Neighbor<Dim<2> >& neighbor = position[nodeListi]->nodeListPtr()->neighbor();
      for (unsigned i = 0; i != n; ++i) {

        const Vector& ri = position(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const Scalar rhoi = rho(nodeListi, i);
        const Vector& gradRhoi = gradRho(nodeListi, i);
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
            phi = min(phi, max(0.0, rij.dot(gradRhoi)*safeInv(rhoi - rhoj)));
          }
        }
        std::sort(pairPlanes.begin(), pairPlanes.end(), compareR2Dplanes);

        // Choose our seed cell shape.
        r2d_poly celli;
        if (numBounds == numNodeLists) {

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

        // Clip the local cell.
        r2d_clip(&celli, &pairPlanes[0], pairPlanes.size());

        // Check if the final polygon is entirely within our "interior" check radius.
        // While we're at it find the average radius (in eta space) of the polygon.
        bool interior = true;
        double Rpoly = 0.0;
        for (unsigned k = 0; k != celli.nverts; ++k) {
          const double Ri = (Hi*Vector(celli.verts[k].pos.x, celli.verts[k].pos.y)).magnitude();
          Rpoly += Ri;
          interior = interior and (Ri < rin);
        }
        Rpoly /= 2.0*celli.nverts;

        if (interior) {

          // This is an interior point -- extract the area.
          r2d_reduce(&celli, voli, 0);
          vol(nodeListi, i) = voli[0];

          // Compute the mass weighted centroid.
          firstmom[0] = rhoi;
          firstmom[1] = phi*gradRhoi.x();
          firstmom[2] = phi*gradRhoi.y();
          r2d_reduce(&celli, firstmom, 1);
          const Scalar m0 = voli[0]*rhoi;
          const Vector deltaCentroidi = Hinv*Vector(firstmom[1], firstmom[2])/m0;

          // If there is a measurable gradient compute the median.
          // Note we convert the gradient here for eta space *and* a unit circle.
          const Vector gradRho_eta = phi*(Hinv*gradRhoi);
          const Scalar b = phi*gradRho_eta.magnitude() * Rpoly;
          Vector deltaMediani;
          if (std::abs(b) >= 0.025*rhoi) {
            const Vector gradUnit = gradRho_eta.unitVector();
            const double tol = 1.0e-5*M_PI*rhoi;
            deltaMediani = Hinv*gradUnit*newtonRaphson(CircleMassAndGradient(rhoi, b), 0.0, 1.0, tol, tol);

            // Combine the centroidal and medial movement.  We take the full medial motion and add the
            // orthogonal bit from the centroid.
            // deltaMedian(nodeListi, i) = 0.5*(deltaMediani + deltaCentroidi);
            deltaMedian(nodeListi, i) = deltaMediani + deltaCentroidi - deltaCentroidi.dot(deltaMediani.unitVector())*deltaCentroidi.unitVector();
            
          } else {

            // Otherwise just use the centroid.
            deltaMedian(nodeListi, i) = deltaCentroidi;

          }

        } else {

          // This point touches a free boundary, so flag it.
          surfacePoint(nodeListi, i) = 1;

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
          facetIndices.back()[1] = 0;
          cells(nodeListi, i) = FacetedVolume(verts, facetIndices);
        }

      }
    }
  }
}

}
}
