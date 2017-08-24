//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/PairComparisons.hh"
#include "Utilities/FastMath.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using namespace FastMath;

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
                     const std::vector<std::vector<Dim<1>::FacetedVolume> >& holes,
                     const FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& weight,
                     const FieldSpace::FieldList<Dim<1>, int>& voidPoint,
                     FieldSpace::FieldList<Dim<1>, int>& surfacePoint,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::Scalar>& vol,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& deltaMedian,
                     FieldSpace::FieldList<Dim<1>, vector<Dim<1>::Vector>>& etaVoidPoints,
                     FieldSpace::FieldList<Dim<1>, Dim<1>::FacetedVolume>& cells) {

  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::SymTensor SymTensor;
  typedef Dim<1>::FacetedVolume FacetedVolume;

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();
  const unsigned numBounds = boundaries.size();
  const bool haveBoundaries = numBounds == numNodeLists;
  const bool haveWeights = weight.size() == numNodeLists;

  REQUIRE(numBounds == 0 or numBounds == numNodeLists);

  const Scalar rin = 0.5*kernelExtent;

  // Zero out return fields.
  surfacePoint = 0;
  deltaMedian = Vector::zero;
  etaVoidPoints = vector<Vector>();

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
  unsigned nodeListj1 = 0, nodeListj2 = 0, j1 = 0, j2 = 0;
  Scalar Hi, H1, H2, rhoi, rho1, rho2, gradRhoi, x1, x2, xi, etamax, b, xm1, xm2, thpt, weighti, weightj, wij,
    xbound0 = -std::numeric_limits<Scalar>::max(),
    xbound1 =  std::numeric_limits<Scalar>::max();

  // cerr << "Volume loop: " << endl;

  // Now walk our sorted point and set the volumes and surface flags.
  
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
      xi = position(nodeListi, i).x();
      Hi = H(nodeListi, i).xx();
      rhoi = rho(nodeListi, i);
      gradRhoi = gradRho(nodeListi, i).x();
      weighti = haveWeights ? weight(nodeListi, i) : 1.0;

      // phi = 1.0;
      if (itr == coords.begin()) {
        x1 = xbound0 - xi;
        H1 = Hi;
        rho1 = rhoi;
        surfacePoint(nodeListi, i) |= 1;
        etaVoidPoints(nodeListi, i).push_back(-rin);
      } else {
        nodeListj1 = (itr-1)->second.first;
        j1 = (itr-1)->second.second;
        H1 = H(nodeListj1, j1).xx();
        rho1 = rho(nodeListj1, j1);
        weightj = haveWeights ? weight(nodeListj1, j1) : 1.0;
        wij = weighti/(weighti + weightj);
        // the direction is here reversed from 2d, so it should weight on the i side
        x1 = wij*(position(nodeListj1, j1).x() - position(nodeListi, i).x());
        // phi = min(phi, max(0.0, gradRhoi*2.0*x1*safeInvVar(rho1 - rhoi)));
        if (voidPoint(nodeListj1, j1) == 1) {
          surfacePoint(nodeListi, i) |= 1;
          etaVoidPoints(nodeListi, i).push_back(-rin);
        } else if (nodeListj1 != nodeListi) {
          surfacePoint(nodeListi, i) |= (1 << (nodeListj1 + 1));
        }
      }

      if (itr == coords.end()-1) {
        x2 = xbound1 - xi;
        H2 = Hi;
        rho2 = rhoi;
        surfacePoint(nodeListi, i) |= 1;
        etaVoidPoints(nodeListi, i).push_back(rin);
      } else {
        nodeListj2 = (itr+1)->second.first;
        j2 = (itr+1)->second.second;
        H2 = H(nodeListj2, j2).xx();
        rho2 = rho(nodeListj2, j2);
        x2 = 0.5*(position(nodeListj2, j2).x() - position(nodeListi, i).x());
        // phi = min(phi, max(0.0, gradRhoi*2.0*x2*safeInvVar(rho2 - rhoi)));
        if (voidPoint(nodeListj2, j2) == 1) {
          surfacePoint(nodeListi, i) |= 1;
          etaVoidPoints(nodeListi, i).push_back(rin);
        } else if (nodeListj2 != nodeListi) {
          surfacePoint(nodeListi, i) |= (1 << (nodeListj2 + 1));
        }
      }

      CHECK(x1 <= 0.0 and x2 >= 0.0);
      // CHECK(phi >= 0.0 and phi <= 1.0);
      etamax = Hi*max(-x1, x2);
      if (surfacePoint(nodeListi, i) == 0 and etamax < rin) {
        vol(nodeListi, i) = x2 - x1;

        b = gradRhoi;
        // cerr << " --> " << i << " " << j1 << " " << j2 << " " << gradRhoi << " " << phi << " " << b << " " << x1 << " " << x2 << endl;
        
        if (std::abs(b)*(x2 - x1) > 1e-8*rhoi) {

          // This version uses the medial position.
          thpt = sqrt(abs(rhoi*rhoi + rhoi*b*(x1 + x2) + 0.5*b*b*(x1*x1 + x2*x2)));
          xm1 = -(rhoi + thpt)/b;
          xm2 = (-rhoi + thpt)/b;
          if (xm1 >= x1 and xm1 <= x2) {
            deltaMedian(nodeListi, i).x(xm1);
          } else {
            deltaMedian(nodeListi, i).x(xm2);
          }
          // cerr << "BLAGO: " << xi << " " << x1 << " " << x2 << " " << b << " " << xm1 << " " << xm2 << " " << deltaMedian(nodeListi, i).x() << " :: "
          //      << rhoi*(xm2 - x1) + 0.5*b*(xm2*xm2 - x1*x1) << " "
          //      << rhoi*(x2 - x1) + 0.5*b*(x2*x2 - x1*x1) << " "
          //      << (rhoi*(xm2 - x1) + 0.5*b*(xm2*xm2 - x1*x1))/(rhoi*(x2 - x1) + 0.5*b*(x2*x2 - x1*x1)) << " "
          //      << endl;

          // This version simply tries rho^2 weighting.
          //deltaMedian(nodeListi, i).x((0.5*rhoi*(x2*x2 - x1*x1) +
          //                             2.0/3.0*rhoi*b*(x2*x2*x2 - x1*x1*x1) +
          //                             0.25*b*b*(x2*x2*x2*x2 - x1*x1*x1*x1))/
          //                            (pow3(rhoi + b*x2) - pow3(rhoi + b*x1)/(3.0*b)));

        } else {
          // Fall back to the straight-up centroid.
          deltaMedian(nodeListi, i).x(0.5*(x2 - x1));
        }

        // Check if the candidate motion is still in the boundary.  If not, project back.
        if (haveBoundaries) {
          const Vector ri = Vector(xi);
          if (not boundaries[nodeListi].contains(ri + deltaMedian(nodeListi, i))) {
            deltaMedian(nodeListi, i) = boundaries[nodeListi].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
          }
          for (unsigned ihole = 0; ihole != holes[nodeListi].size(); ++ihole) {
            if (holes[nodeListi][ihole].contains(ri + deltaMedian(nodeListi, i))) {
              deltaMedian(nodeListi, i) = holes[nodeListi][ihole].closestPoint(ri + deltaMedian(nodeListi, i)) - ri;
            }
          }
        }

        // {
        //   const Scalar x0 = -dx1,
        //     x1 = dx2,
        //     xm = deltaMedian(nodeListi, i).x(),
        //     m1 = -rhoi*x0 + 0.5*b*x0*x0 + rhoi*xm + 0.5*b*xm*xm,
        //     m2 = -rhoi*xm + 0.5*b*xm*xm + rhoi*x1 + 0.5*b*x1*x1;
        //   cerr << " --> " << i << " " << x0 << " " << x1 << " " << xm << " " << m1 << " " << m2 << endl;
        // }

      } else {
        surfacePoint(nodeListi, i) |= 1;
        if (-Hi*x1 >= rin) etaVoidPoints(nodeListi, i).push_back(max(Hi*x1, -rin));
        if ( Hi*x2 >= rin) etaVoidPoints(nodeListi, i).push_back(min(Hi*x2,  rin));
      }
    }
    CHECK2(((surfacePoint(nodeListi, i) & 1) == 1 and
            (etaVoidPoints(nodeListi, i).size() == 1 or etaVoidPoints(nodeListi, i).size() == 2)) or
           (surfacePoint(nodeListi, i) & 1) == 0,
           "(" << nodeListi << " " << i << ") " << surfacePoint(nodeListi, i) << " " << etaVoidPoints(nodeListi, i).size());

    // cerr << "  " << i << " " << vol(nodeListi, i) << " " << surfacePoint(nodeListi, i) << " "
    //      << j1 << " " << j2 << " " << voidPoint(nodeListj1, j1) << " " << voidPoint(nodeListj2, j2)
    //      << " ---- " << position(nodeListj1, j1).x() << " " << position(nodeListi, i) << " " << position(nodeListj2, j2).x() 
    //      << endl;

  }

}

}
}
