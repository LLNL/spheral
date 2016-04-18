//---------------------------------Spheral++------------------------------------
// Compute the volume per point based on the Voronoi tessellation.
//------------------------------------------------------------------------------
#include "computeVoronoiVolume.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/allReduce.hh"
#include "Mesh/Mesh.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;

using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
void
computeVoronoiVolume(const FieldList<Dim<1>, Dim<1>::Vector>& position,
                     FieldList<Dim<1>, Dim<1>::Scalar>& vol) {

  const unsigned numGens = position.numNodes();
  const unsigned numNodeLists = position.size();

  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::SymTensor SymTensor;
  typedef Dim<1>::FacetedVolume FacetedVolume;

  // Copy the input positions to single list.
  // For convenience we gather all ghost node coordinates at the end of the
  // list.
  Vector xmin(1e100), xmax(-1e100);
  vector<Vector> coords;
  coords.reserve(numGens);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = position[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      coords.push_back(position(nodeListi, i));
      const Scalar xi = position(nodeListi, i).x();
      xmin.x(min(xmin.x(), xi));
      xmax.x(max(xmax.x(), xi));
    }
  }
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned firstGhostNode = position[nodeListi]->nodeListPtr()->firstGhostNode();
    const unsigned n = position[nodeListi]->nodeListPtr()->numGhostNodes();
    for (unsigned i = 0; i != n; ++i) {
      coords.push_back(position(nodeListi, firstGhostNode + i));
      const Scalar xi = position(nodeListi, firstGhostNode + i).x();
      xmin.x(min(xmin.x(), xi));
      xmax.x(max(xmax.x(), xi));
    }
  }
  xmin.x(allReduce(xmin.x(), MPI_MIN, Communicator::communicator()));
  xmax.x(allReduce(xmax.x(), MPI_MAX, Communicator::communicator()));
  const Vector delta = 0.01*(xmax - xmin);
  xmin -= delta;
  xmax += delta;
  CHECK(coords.size() == numGens);

  // Do the tessellation.
  const MeshSpace::Mesh<Dim<1> > mesh(coords, xmin, xmax);

  // Now we can extract the volumes.
  unsigned icell = 0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = vol[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i, ++icell) {
      vol(nodeListi, i) = mesh.zone(icell).volume();
    }
  }
}

}
}
