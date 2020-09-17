//---------------------------------Spheral++----------------------------------//
// copy2polytope
//
// Helper method to copy a set of Spheral polyhedra to a polytope tessellation.
//
// Created by JMO, Wed Jan  9 16:25:28 PST 2019
//----------------------------------------------------------------------------//
#include "copy2polytope.hh"

namespace Spheral {

template<typename Dimension, int nDim>
void copy2polytope(const FieldList<Dimension, typename Dimension::FacetedVolume>& cells,
                   polytope::Tessellation<nDim, double>& mesh) {
                   
  mesh.clear();

  const auto numNodeLists = cells.size();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = cells[nodeListi]->numInternalElements();
    const auto noldcells = mesh.cells.size();
    mesh.cells.resize(noldcells + n);
    for (auto i = 0u; i < n; ++i) {
      const auto& celli = cells(nodeListi, i);
      const auto& verts = celli.vertices();
      const auto& facets = celli.facets();
      const auto  noldnodes = mesh.nodes.size()/nDim;
      const auto  noldfaces = mesh.faces.size();
      mesh.faces.resize(noldfaces + facets.size());
      for (auto j = 0u; j < verts.size(); ++j) {
        for (auto k = 0; k < nDim; ++k) mesh.nodes.push_back(verts[j][k]);
      }
      for (auto j = 0u; j < facets.size(); ++j) {
        mesh.cells[noldcells + i].push_back(noldfaces + j);
        const auto& ipoints = facets[j].ipoints();
        for (auto k: ipoints) {
          mesh.faces[noldfaces + j].push_back(noldnodes + k);
        }
      }
    }
  }
}

}
