#ifndef _Spheral_NeighbourSpace_NodePairList_hh_
#define _Spheral_NeighbourSpace_NodePairList_hh_

#include <iostream>
#include <vector>

namespace Spheral {

struct NodePairIdxType {
  NodePairIdxType(int i_n, int i_l, int j_n, int j_l);
  int i_node, i_list, j_node, j_list;
};


class NodePairList {
public:
  typedef std::vector<NodePairIdxType> ContainerType;
  typedef typename ContainerType::value_type value_type;
  typedef typename ContainerType::reference reference;
  typedef typename ContainerType::const_reference const_reference;
  typedef typename ContainerType::iterator iterator;
  typedef typename ContainerType::const_iterator const_iterator;
  typedef typename ContainerType::reverse_iterator reverse_iterator;
  typedef typename ContainerType::const_reverse_iterator const_reverse_iterator;

  NodePairList();
  void push_back(NodePairIdxType nodePair);
  void clear(); 
  size_t size() const { return mNodePairList.size(); }

  iterator begin() { return mNodePairList.begin(); }
  iterator end() { return mNodePairList.end(); }
  const_iterator begin() const { return mNodePairList.begin(); }
  const_iterator end() const { return mNodePairList.end(); }

  reverse_iterator rbegin() { return mNodePairList.rbegin(); }
  reverse_iterator rend() { return mNodePairList.rend(); }
  const_reverse_iterator rbegin() const { return mNodePairList.rbegin(); }
  const_reverse_iterator rend() const { return mNodePairList.rend(); }

  reference operator[](const size_t i) { return mNodePairList[i]; }
  const_reference operator[](const size_t i) const { return mNodePairList[i]; }

  template<typename InputIterator>
  iterator insert(const_iterator pos, InputIterator first, InputIterator last) { return mNodePairList.insert(pos, first, last); }

private:
  ContainerType mNodePairList;
};


} //namespace Spheral


#endif // _Spheral_NeighbourSpace_NodePairList_hh_
