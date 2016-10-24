//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/PairComparisons.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
void
computeVoronoiVolume(const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& position,
                     const FieldSpace::FieldList<Dim<1>, Dim<1>::SymTensor>& H,
                     const FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& rho,
                     const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& gradRho,
                     const NeighborSpace::ConnectivityMap<Dim<1> >& connectivityMap,
                     const Dim<1>::Scalar kernelExtent,
                     FieldSpace::FieldList<Dim<1>, int>& surfacePoint,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& deltaCentroid) {

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();

  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::SymTensor SymTensor;
  typedef Dim<1>::FacetedVolume FacetedVolume;

  const Scalar rin = 0.5*kernelExtent;

  // Zero out the deltaCentroid field.
  deltaCentroid = Vector::zero;

  // Copy the input positions to single list, and sort it.
  // Note our logic here relies on ghost nodes already being built, including parallel nodes.
  typedef pair<double, pair<unsigned, unsigned> > PointCoord;
  vector<PointCoord> coords;
  coords.reserve(numGens);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = position[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      coords.push_back(make_pair(position(nodeListi, i).x(), make_pair(nodeListi, i)));
    }
  }
  sort(coords.begin(), coords.end(), ComparePairsByFirstElement<PointCoord>());

  // Now walk our sorted point and set the volumes and surface flags.
  surfacePoint = 0;
  const vector<NodeList<Dim<1> >*>& nodeListPtrs = position.nodeListPtrs();
  for (vector<PointCoord>::const_iterator itr = coords.begin();
       itr != coords.end();
       ++itr) {
    const unsigned nodeListi = itr->second.first;
    const unsigned i = itr->second.second;
    if (i < nodeListPtrs[nodeListi]->firstGhostNode()) {
      if (itr == coords.begin() or itr == coords.end()-1) {
        surfacePoint(nodeListi, i) = 1;
      } else {
        const unsigned nodeListj1 = (itr-1)->second.first,
                       nodeListj2 = (itr+1)->second.first,
                               j1 = (itr-1)->second.second,
                               j2 = (itr+1)->second.second;
        const Scalar Hi = H(nodeListi, i).xx(),
                     H1 = H(nodeListj1, j1).xx(),
                     H2 = H(nodeListj2, j2).xx(),
                   rhoi = rho(nodeListi, i),
                   rho1 = rho(nodeListj1, j1),
                   rho2 = rho(nodeListj2, j2),
               gradRhoi = gradRho(nodeListi, i).x();
        const Scalar dx1 = position(nodeListi, i).x() - position(nodeListj1, j1).x(),
                     dx2 = position(nodeListj2, j2).x() - position(nodeListi, i).x();
        CHECK(dx1 >= 0.0 and dx2 >= 0.0);
        const Scalar etamin = min(Hi, min(H1, H2))*min(dx1, dx2);
        if (etamin < rin) {
          vol(nodeListi, i) = 0.5*(dx1 + dx2);
          const Scalar phi = min(1.0, min(max(0.0, dx1*safeInvVar(rhoi - rho1)),
                                          max(0.0, dx2*safeInvVar(rho2 - rhoi))));
          CHECK(phi >= 0.0 and phi <= 1.0);
          const Scalar m1 = 0.5*dx1*(rhoi - phi*gradRhoi*0.25*dx1),
                       m2 = 0.5*dx2*(rhoi + phi*gradRhoi*0.25*dx2);
          deltaCentroid(nodeListi, i).x(0.5*(m2*dx2 - m1*dx1)/(m1 + m2));
        } else {
          surfacePoint(nodeListi, i) = 1;
        }
      }
    }
  }
}

}
}
