#ifndef _Spheral_NeighbourSpace_NodePairList_hh_
#define _Spheral_NeighbourSpace_NodePairList_hh_

#include <iostream>


namespace Spheral {

class NodePairIdxType {
  public:
    NodePairIdxType(int in, int il, int jn, int jl) : 
      i_node(in), i_list(il), j_node(jn), j_list(jl) {}
  private:
    int i_node, i_list, j_node, j_list;
};


class NodePairList {
  public:
    NodePairList(){ }
    void push_back(NodePairIdxType nodePair) { mNodePairList.push_back(nodePair); }
    void clear() { mNodePairList.clear(); }
    unsigned int size() const { return (unsigned int)mNodePairList.size(); }
  private:
    std::vector<NodePairIdxType> mNodePairList;
};


} //namespace Spheral


#endif // _Spheral_NeighbourSpace_NodePairList_hh_
