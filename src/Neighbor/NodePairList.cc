#include "NodePairList.hh"

namespace Spheral {
  
  NodePairIdxType::NodePairIdxType(int i_n, int i_l, int j_n, int j_l) :
    i_node(i_n), i_list(i_l), j_node(j_n), j_list(j_l) {}

  NodePairList::NodePairList(){};

  void NodePairList::push_back(NodePairIdxType nodePair) {
    mNodePairList.push_back(nodePair);
  } 

  void NodePairList::clear() {
    mNodePairList.clear();
  }

}
