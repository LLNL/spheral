//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
extern "C" {
#include "r3d/r2d.h"
}

#include <algorithm>

#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/allReduce.hh"

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
                     FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& deltaMedian) {

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();
  const unsigned numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());

  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  if (numGensGlobal > 0) {

    // Start out assuming all points are internal.
    surfacePoint = 0;

    // Build an approximation of the starting kernel shape.
    const unsigned nverts = 18;
    const Scalar rin2 = 0.25*kernelExtent*kernelExtent * FastMath::square(cos(M_PI/nverts));
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

            // Build the half plane.
            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector nhat = etai.unitVector();
            pairPlanes.push_back(r2d_plane());
            pairPlanes.back().n.x = nhat.x();
            pairPlanes.back().n.y = nhat.y();
            pairPlanes.back().d = 0.5*etai.magnitude();

            // Check the density gradient limiter.
            phi = min(phi, max(0.0, rij.dot(gradRhoi)*safeInv(rhoi - rhoj)));
          }
        }
        std::sort(pairPlanes.begin(), pairPlanes.end(), compareR2Dplanes);

        // Start with the initial cell shape (in eta space).
        r2d_poly celli = initialCell;
        CHECK2(r2d_is_good(&celli), "Bad polygon!");

        // Clip the local cell.
        r2d_clip(&celli, &pairPlanes[0], pairPlanes.size());

        // Are there any of the original volume vertices left?
        bool interior = true;
        {
          unsigned k = 0;
          do {
            interior = (FastMath::square(celli.verts[k].pos.x) + FastMath::square(celli.verts[k].pos.y) < rin2);
          } while (interior and ++k != celli.nverts);
        }
        if (interior) {

          // This is an interior point -- extract the area.
          r2d_reduce(&celli, voli, 0);
          vol(nodeListi, i) = voli[0]/Hdeti;

          // Also get the centroid.
          // Note we have to convert the gradient to eta space, since that's what the polygon is defined in.
          const SymTensor Hinv = Hi.Inverse();
          const Vector gradRho_eta = phi*(Hinv*gradRhoi);
          firstmom[0] = rhoi;
          firstmom[1] = gradRho_eta.x();
          firstmom[2] = gradRho_eta.y();
          r2d_reduce(&celli, firstmom, 1);
          const Scalar m0 = voli[0]*rhoi;
          deltaMedian(nodeListi, i) = Hinv*Vector(firstmom[1], firstmom[2])/m0;

        } else {

          // This point touches a free boundary, so flag it.
          surfacePoint(nodeListi, i) = 1;

        }
      }
    }
  }
}

}
}
