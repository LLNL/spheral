#ifndef _Spheral_NeighbourSpace_NodePairList_hh_
#define _Spheral_NeighbourSpace_NodePairList_hh_

#include <iostream>


namespace Spheral {

class NodePairIdxType {
  public:
    NodePairIdxType(int i_n, int i_l, int j_n, int j_l);
  private:
    int i_node, i_list, j_node, j_list;
};


class NodePairList {
  public:
    NodePairList();
    void push_back(NodePairIdxType nodePair);
    void clear(); 
    unsigned int size() const; 
  private:
    std::vector<NodePairIdxType> mNodePairList;
};


} //namespace Spheral


#endif // _Spheral_NeighbourSpace_NodePairList_hh_
