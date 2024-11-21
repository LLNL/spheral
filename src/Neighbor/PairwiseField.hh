//---------------------------------Spheral++----------------------------------//
// PairwiseField
//
// Stores a value per node pair in a NodePairList. Because connectivity is
// allowed to change step to step in our meshfree methods, PairwiseField is
// ephemeral and will be invalidated when topology is updated.
//
// Created by J. Michael Owen, Wed Nov 20 14:44:44 PST 2024
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_PairwiseField_hh_
#define _Spheral_NeighborSpace_PairwiseField_hh_

#include "Utilities/DBC.hh"

#include <vector>
#include <memory>

namespace Spheral {

// Forward declarations
template<typename Dimension> class ConnectivityMap;
class NodePairList;

template<typename Dimension, typename Value>
class PairwiseField {
public:
  //--------------------------- Public Interface ---------------------------//
  using ContainerType = std::vector<Value>;
  using value_type = typename ContainerType::value_type;
  using reference = typename ContainerType::reference;
  using const_reference = typename ContainerType::const_reference;
  using iterator = typename ContainerType::iterator;
  using const_iterator = typename ContainerType::const_iterator;
  using reverse_iterator = typename ContainerType::reverse_iterator;
  using const_reverse_iterator = typename ContainerType::const_reverse_iterator;

  // Constructors, destructors
  PairwiseField(const ConnectivityMap<Dimension>& connectivity);
  PairwiseField(const PairwiseField& rhs)                                  = default;
  ~PairwiseField()                                                         = default;
  PairwiseField& operator=(const PairwiseField& rhs)                       = default;

  // Access the data
  const Value& operator[](const size_t k) const                            { REQUIRE(k < mValues.size()); return mValues[k]; }
  const Value& operator()(const size_t k) const                            { return (*this)[k]; }
  const Value& operator()(const NodePairIdxType& x) const;

  Value& operator[](const size_t k)                                        { REQUIRE(k < mValues.size()); return mValues[k]; }
  Value& operator()(const size_t k)                                        { return (*this)[k]; }
  Value& operator()(const NodePairIdxType& x);

  // Iterators
  const_iterator begin() const                                             { return mValues.begin(); }
  const_iterator end() const                                               { return mValues.end(); }
  const_iterator rbegin() const                                            { return mValues.rbegin(); }
  const_iterator rend() const                                              { return mValues.rend(); }

  iterator begin()                                                         { return mValues.begin(); }
  iterator end()                                                           { return mValues.end(); }
  iterator rbegin()                                                        { return mValues.rbegin(); }
  iterator rend()                                                          { return mValues.rend(); }

  // Other methods
  const NodePairList& pairs() const;
  size_t size() const                                                      { return mValues.size(); }

  // Forbidden methods
  PairwiseField()                                                          = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  std::weak_ptr<NodePairList> mPairsPtr;
  ContainerType mValues;
};

}

#include "Neighbor/PairwiseFieldInline.hh"

#endif

