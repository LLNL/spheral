#ifndef DomainNode_HH
#define DomainNode_HH

#include <stddef.h>

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

namespace Spheral {
namespace PartitionSpace {

template<typename Dimension>
struct DomainNode {
  int localNodeID;
  int uniqueLocalNodeID;
  int globalNodeID;
  int nodeListID;
  int domainID;
  double work;
  typename Dimension::Vector position;

  // These methods are required by the Boost.Python vector_indexing_suite.
  bool operator==(const DomainNode& rhs) const;
  bool operator!=(const DomainNode& rhs) const;

  // Helpful methods for parallel communication of DomainNodes.
  static size_t packSize();
  std::vector<double> pack() const;
  void unpack(std::vector<double>::const_iterator& itr);
};

}
}

#ifndef __GCCXML__
#include "DomainNodeInline.hh"
#endif

#else
namespace Spheral {
  namespace PartitionSpace {
    template<typename Dimension> struct DomainNode;
  }
}

#endif
