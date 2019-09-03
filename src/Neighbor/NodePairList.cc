#include "NodePairList.hh"

namespace Spheral {
  
  NodePairIdxType::NodePairIdxType(int i_n, int i_l, int j_n, int j_l) :
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
