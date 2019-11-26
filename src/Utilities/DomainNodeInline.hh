#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// operator ==
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
DomainNode<Dimension>::
operator==(const DomainNode<Dimension>& rhs) const {
  return (localNodeID == rhs.localNodeID &&
          uniqueLocalNodeID == rhs.uniqueLocalNodeID &&
          globalNodeID == rhs.globalNodeID &&
          nodeListID == rhs.nodeListID &&
          domainID == rhs.domainID &&
          position == rhs.position &&
          work == rhs.work);
}

//------------------------------------------------------------------------------
// operator !=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
DomainNode<Dimension>::
operator!=(const DomainNode<Dimension>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// The size of the packed representation of a DomainNode.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
size_t
DomainNode<Dimension>::
packSize() {
  return 6 + Dimension::nDim;
}

//------------------------------------------------------------------------------
// Pack the domain node into a std::vector<double> representation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<double>
DomainNode<Dimension>::
pack() const {
  std::vector<double> result;
  result.reserve(packSize());
  result.push_back(double(localNodeID));
  result.push_back(double(uniqueLocalNodeID));
  result.push_back(double(globalNodeID));
  result.push_back(double(nodeListID));
  result.push_back(double(domainID));
  result.push_back(double(work));
  for (typename Dimension::Vector::const_iterator itr = position.begin();
       itr != position.end();
       ++itr) result.push_back(*itr);
  ENSURE(result.size() == packSize());
  return result;
}

//------------------------------------------------------------------------------
// Unpack the domain nodes state from a std::vector<double> representation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
DomainNode<Dimension>::
unpack(std::vector<double>::const_iterator& itr) {
  localNodeID = int(*(itr++));
  uniqueLocalNodeID = int(*(itr++));
  globalNodeID = int(*(itr++));
  nodeListID = int(*(itr++));
  domainID = int(*(itr++));
  work = *(itr++);
  for (int i = 0; i != Dimension::nDim; ++i) position(i) = *(itr++);
}

}
