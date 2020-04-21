#ifndef _Spheral_NeighbourSpace_NodePairList_hh_
#define _Spheral_NeighbourSpace_NodePairList_hh_

#include <iostream>
#include <vector>

namespace Spheral {

struct NodePairIdxType {
  NodePairIdxType(int i_n, int i_l, int j_n, int j_l,
                  double f = 1.0);
  int i_node, i_list, j_node, j_list;
  double f_couple;                       // An arbitrary fraction in [0,1] to hold the effective coupling of the pair

  // Comparisons
  bool operator==(const NodePairIdxType& val) const { return (i_list == val.i_list and
                                                              i_node == val.i_node and
                                                              j_list == val.j_list and
                                                              j_node == val.j_node); }
  bool operator<(const NodePairIdxType& val) const {  return (i_list < val.i_list ? true :
                                                              i_list > val.i_list ? false :
                                                              i_node < val.i_node ? true :
                                                              i_node > val.i_node ? false :
                                                              j_list < val.j_list ? true :
                                                              j_list > val.j_list ? false :
                                                              j_node < val.j_node ? true :
                                                              false); }

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

  // Iterators
  iterator begin() { return mNodePairList.begin(); }
  iterator end() { return mNodePairList.end(); }
  const_iterator begin() const { return mNodePairList.begin(); }
  const_iterator end() const { return mNodePairList.end(); }

  // Reverse iterators
  reverse_iterator rbegin() { return mNodePairList.rbegin(); }
  reverse_iterator rend() { return mNodePairList.rend(); }
  const_reverse_iterator rbegin() const { return mNodePairList.rbegin(); }
  const_reverse_iterator rend() const { return mNodePairList.rend(); }

  // Indexing
  reference operator[](const size_t i) { return mNodePairList[i]; }
  const_reference operator[](const size_t i) const { return mNodePairList[i]; }

  // Inserting
  template<typename InputIterator>
  iterator insert(const_iterator pos, InputIterator first, InputIterator last) { return mNodePairList.insert(pos, first, last); }

private:
  ContainerType mNodePairList;
};


} //namespace Spheral


#endif // _Spheral_NeighbourSpace_NodePairList_hh_
