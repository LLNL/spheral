//------------------------------------------------------------------------------
// generateVoidNodes
//
// This algorithm tries to analyze how continuous a node distribution is, and 
// if it determines there is an edge to the distribution creates new void nodes
// outside that surface.
// We assume here that the caller has already created all the boundary ghost 
// nodes.
//------------------------------------------------------------------------------
#ifndef __Spheral_generateVoidNodes__
#define __Spheral_generateVoidNodes__

#include <vector>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class DataBase;
template<typename Dimension> class Mesh;
template<typename Dimension> class NodeList;

// Our method.
template<typename Dimension>
void generateVoidNodes(const std::vector<typename Dimension::Vector>& generators,
                       const std::vector<typename Dimension::SymTensor>& Hs,
                       const Mesh<Dimension>& mesh,
                       const typename Dimension::Vector& xmin,
                       const typename Dimension::Vector& xmax,
                       const unsigned numInternal,
                       const double nPerh,
                       const double threshold,
                       NodeList<Dimension>& voidNodes);

}

#endif
