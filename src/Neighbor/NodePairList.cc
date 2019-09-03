#include "NodePairList.hh"

namespace Spheral {
  
  NodePairIdxType::NodePairIdxType(int in, int il, int jn, int jl) :
    i_node(in), i_list(il), j_node(jn), j_list(jl) {}

  NodePairList::NodePairList(){};

  void NodePairList::push_back(NodePairIdxType nodePair) {
    mNodePairList.push_back(nodePair);
  } 

  void NodePairList::clear() {
    mNodePairList.clear();
  }

  unsigned int NodePairList::size() const {
    return (unsigned int)mNodePairList.size();
  }

}
