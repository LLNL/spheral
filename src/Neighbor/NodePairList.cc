#include "NodePairList.hh"

#include <algorithm>

namespace Spheral {
  
  NodePairIdxType::NodePairIdxType(int i_n, int i_l, int j_n, int j_l, double f) :
    i_node(i_n), i_list(i_l), j_node(j_n), j_list(j_l), f_couple(f) {}

  NodePairList::NodePairList(){}

  void NodePairList::push_back(NodePairIdxType nodePair) {
    mNodePairList.push_back(nodePair);
  } 

  void NodePairList::clear() {
    mNodePairList.clear();
  }

  void NodePairList::makeUnique() {
    // Make sure the node pairs are ordered correctly
    for (auto kk = 0u; kk < mNodePairList.size(); ++kk)
    {
      auto& pair = mNodePairList[kk];
      if (pair.i_list > pair.j_list)
      {
        const auto temp_list = pair.j_list;
        const auto temp_node = pair.j_node;
        pair.j_list = pair.i_list;
        pair.j_node = pair.i_node;
        pair.i_list = temp_list;
        pair.i_node = temp_node;
      }
      if (pair.i_list == pair.j_list && pair.i_node > pair.j_node)
      {
        const auto temp = pair.j_node;
        pair.j_node = pair.i_node;
        pair.i_node = temp;
      }
    }
    
    // Remove duplicates
    std::sort(mNodePairList.begin(), mNodePairList.end());
    mNodePairList.erase(std::unique(mNodePairList.begin(), mNodePairList.end()), mNodePairList.end());
  }

}
