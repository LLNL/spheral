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
                     const std::vector<Dim<1>::FacetedVolume>& boundaries,
                     FieldSpace::FieldList<Dim<1>, int>& surfacePoint,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& deltaMedian,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells) {

  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::SymTensor SymTensor;
  typedef Dim<1>::FacetedVolume FacetedVolume;

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();
  const unsigned numBounds = boundaries.size();

  REQUIRE(numBounds == 0 or numBounds == numNodeLists);

  const Scalar rin = 0.5*kernelExtent;

  // Zero out the deltaMedian field.
  deltaMedian = Vector::zero;

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

  // Prepare some scratch variables.
  unsigned nodeListj1, nodeListj2, j1, j2;
  Scalar Hi, H1, H2, rhoi, rho1, rho2, gradRhoi, dx1, dx2, etamax, phi, b, rho0, dx, m1, m2,
    xbound0 = -std::numeric_limits<Scalar>::max(),
    xbound1 =  std::numeric_limits<Scalar>::max();

  // Now walk our sorted point and set the volumes and surface flags.
  surfacePoint = 0;
  const vector<NodeList<Dim<1> >*>& nodeListPtrs = position.nodeListPtrs();
  for (vector<PointCoord>::const_iterator itr = coords.begin();
       itr != coords.end();
       ++itr) {
    const unsigned nodeListi = itr->second.first;
    const unsigned i = itr->second.second;
    if (i < nodeListPtrs[nodeListi]->firstGhostNode()) {

      // Is there a bounding volume for this NodeList?
      if (numBounds > 0) {
        xbound0 = boundaries[nodeListi].xmin().x();
        xbound1 = boundaries[nodeListi].xmax().x();
      }

      // Grab our state.
      Hi = H(nodeListi, i).xx();
      rhoi = rho(nodeListi, i);
      gradRhoi = gradRho(nodeListi, i).x();

      if (itr == coords.begin()) {
        dx1 = position(nodeListi, i).x() - xbound0;
        H1 = Hi;
        rho1 = rhoi; // - gradRhoi*dx1;
      } else {
        nodeListj1 = (itr-1)->second.first;
        j1 = (itr-1)->second.second;
        H1 = H(nodeListj1, j1).xx();
        rho1 = rho(nodeListj1, j1);
        dx1 = position(nodeListi, i).x() - position(nodeListj1, j1).x();
      }

      if (itr == coords.end()-1) {
        dx2 = xbound1 - position(nodeListi, i).x();
        H2 = Hi;
        rho2 = rhoi; //  + gradRhoi*dx2;
      } else {
        nodeListj2 = (itr+1)->second.first;
        j2 = (itr+1)->second.second;
        H2 = H(nodeListj2, j2).xx();
        rho2 = rho(nodeListj2, j2);
        dx2 = position(nodeListj2, j2).x() - position(nodeListi, i).x();
      }

      CHECK(dx1 >= 0.0 and dx2 >= 0.0);
      etamax = max(Hi, max(H1, H2))*max(dx1, dx2);
      if (etamax < rin) {
        vol(nodeListi, i) = 0.5*(dx1 + dx2);
        const Scalar phi = min(1.0, min(max(0.0, dx1*safeInvVar(rhoi - rho1)),
                                        max(0.0, dx2*safeInvVar(rho2 - rhoi))));
        CHECK(phi >= 0.0 and phi <= 1.0);

        const Scalar b = phi*gradRhoi;
        const Scalar rho0 = rhoi - b*dx1;
        const Scalar dx = dx1 + dx2;
        if (std::abs(b*dx) > 0.05*rho0) {
          deltaMedian(nodeListi, i).x((sqrt(2.0*rho0*rho0 + b*b*dx*dx + 2.0*b*rho0*dx)/sqrt(2.0) - rho0)*safeInvVar(b) - dx1);
        } else {
          m1 = 0.5*dx1*(rhoi - phi*gradRhoi*0.25*dx1),
          m2 = 0.5*dx2*(rhoi + phi*gradRhoi*0.25*dx2);
          deltaMedian(nodeListi, i).x(0.5*(m2*dx2 - m1*dx1)/(m1 + m2));
        }

        // {
        //   const Scalar x0 = -0.5*dx1,
        //     x1 = 0.5*dx2,
        //     xm = deltaMedian(nodeListi, i).x(),
        //     m1 = (rhoi + b*0.5*x0)*abs(x0),
        //     m2 = (rhoi + b*0.5*x1)*x1,
        //     mm1 = (rhoi + b*0.5*(x0 + xm))*(xm - x0),
        //     mm2 = (rhoi + b*0.5*(x1 + xm))*(x1 - xm);
        //   cerr << " --> " << i << " " << x0 << " " << x1 << " " << xm << " " << m1 << " " << m2 << " " << mm1 << " " << mm2 << endl;
        // }

      } else {
        surfacePoint(nodeListi, i) = 1;
      }
    }
  }
}

}
}
