//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include "polytope/polytope.hh"

#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/allReduce.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
void
computeVoronoiVolume(const FieldList<Dim<2>, Dim<2>::Vector>& position,
                     FieldList<Dim<2>, Dim<2>::Scalar>& vol) {

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();
  const unsigned numGensGlobal = allReduce(numGens, MPI_SUM, Communicator::communicator());

  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::FacetedVolume FacetedVolume;

  if (numGensGlobal > 0) {

    // Copy the input positions to polytope's convention.
    // For convenience we gather all ghost node coordinates at the end of the
    // list.
    vector<double> coords;
    coords.reserve(2*numGens);
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        const Vector& ri = position(nodeListi, i);
        std::copy(ri.begin(), ri.end(), std::back_inserter(coords));
      }
    }
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned firstGhostNode = position[nodeListi]->nodeListPtr()->firstGhostNode();
      const unsigned n = position[nodeListi]->nodeListPtr()->numGhostNodes();
      for (unsigned i = 0; i != n; ++i) {
        const Vector& ri = position(nodeListi, firstGhostNode + i);
        std::copy(ri.begin(), ri.end(), std::back_inserter(coords));
      }
    }
    CHECK(coords.size() == 2*numGens);

    // Do the polytope tessellation.
    polytope::Tessellation<2, double> tessellation;
    {
#if defined USE_TRIANGLE && ( USE_TRIANGLE==1 )
      polytope::TriangleTessellator<double> tessellator;
#else
      polytope::BoostTessellator<double> tessellator;
#endif
      tessellator.tessellateDegenerate(coords, 1.0e-8, tessellation);
    }
    CHECK(tessellation.cells.size() == numGens);

    // Now we can extract the areas.
    const vector<double>& nodes = tessellation.nodes;
    const vector<vector<unsigned> >& faces = tessellation.faces;
    unsigned icell = 0;
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = vol[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i, ++icell) {
        const Vector& ri = position(nodeListi, i);
        const vector<int>& faceIDs = tessellation.cells[icell];
        const unsigned nfaces = faceIDs.size();
        CHECK(nfaces >= 3);
        Scalar areasum = 0.0;
        for (unsigned j = 0; j != nfaces; ++j) {
          int iface = faceIDs[j] < 0 ? ~faceIDs[j] : faceIDs[j];
          CHECK(faces[iface].size() == 2);
          const Vector a = Vector(nodes[2*faces[iface][0]], nodes[2*faces[iface][0] + 1]);
          const Vector b = Vector(nodes[2*faces[iface][1]], nodes[2*faces[iface][1] + 1]);
          areasum += abs(((a - ri).cross(b - ri)).z());
        }
        vol(nodeListi, i) = 0.5*areasum;
      }
    }
  }
}

}
}
