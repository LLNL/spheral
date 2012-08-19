//---------------------------------Spheral++----------------------------------//
// MortonOrderRedistributeNodes
// Attempt to redistribute nodes such that they are laid out in memory
// in a Morton ordering.  Note that this involves renumbering the nodes of 
// each NodeList, not just redistributing them between processors.
//
// Created by JMO, Tue Mar 25 14:19:18 PDT 2008
//----------------------------------------------------------------------------//
#include "PeanoHilbertOrderRedistributeNodes.hh"
#include "Utilities/peanoHilbertOrderIndicies.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"

namespace Spheral {
namespace PartitionSpace {

using DataBaseSpace::DataBase;
using FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PeanoHilbertOrderRedistributeNodes<Dimension>::
PeanoHilbertOrderRedistributeNodes(double dummy,
                                   const double minNodesPerDomainFraction,
                                   const double maxNodesPerDomainFraction,
                                   const bool workBalance,
                                   const bool localReorderOnly):
  SpaceFillingCurveRedistributeNodes<Dimension>(dummy, 
                                                minNodesPerDomainFraction,
                                                maxNodesPerDomainFraction,
                                                workBalance,
                                                localReorderOnly) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PeanoHilbertOrderRedistributeNodes<Dimension>::
~PeanoHilbertOrderRedistributeNodes() {
}

//------------------------------------------------------------------------------
// Hash the node positions into their tree ordered indicies.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, KeyTraits::Key>
PeanoHilbertOrderRedistributeNodes<Dimension>::
computeHashedIndicies(const DataBase<Dimension>& dataBase) const {
  return peanoHilbertOrderIndicies(dataBase);
}

}
}

