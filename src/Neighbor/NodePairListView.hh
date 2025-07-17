#ifndef Spheral_NodePairListView_hh
#define Spheral_NodePairListView_hh

#include "chai/ManagedArray.hpp"
#include "Neighbor/NodePairIdxType.hh"

namespace Spheral {

//------------------------------------------------------------------------------
class NodePairListView {

  using ContainerType = typename chai::ManagedArray<NodePairIdxType>;
  ContainerType mData;

public:
  SPHERAL_HOST_DEVICE NodePairListView() = default;
  SPHERAL_HOST NodePairListView(ContainerType const &d) : mData(d) {}
  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) { return mData[i]; }
  SPHERAL_HOST_DEVICE
  NodePairIdxType& operator[](const size_t i) const { return mData[i]; }
  SPHERAL_HOST_DEVICE
  size_t size() const { return mData.size(); }
  void move(chai::ExecutionSpace space) { mData.move(space); }
};
}
#endif
