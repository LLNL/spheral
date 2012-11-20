//---------------------------------Spheral++----------------------------------//
// computeGenerators.
//
// Helper method used when we're generating meshes.
// This method takes a set of NodeLists, figures out what domains have to talk
// to one another, and returns the flattened set of positions and Hs for the 
// the generators this domain needs (including those from neighbor domains).
//
// Created by JMO, Mon Dec  6 10:34:43 PST 2010
//----------------------------------------------------------------------------//
#include <vector>

namespace Spheral {
  namespace MeshSpace {
    template<typename Dimension, typename NodeListIterator, typename BoundaryIterator>
    void
    computeGenerators(NodeListIterator nodeListBegin,
                      NodeListIterator nodeListEnd,
                      BoundaryIterator boundaryBegin,
                      BoundaryIterator boundaryEnd,
                      const typename Dimension::Vector& xmin,
                      const typename Dimension::Vector& xmax,
                      const bool generateParallelRind,
                      std::vector<typename Dimension::Vector>& positions,
                      std::vector<typename Dimension::SymTensor>& Hs,
                      std::vector<unsigned>& offsets);
  }
}
