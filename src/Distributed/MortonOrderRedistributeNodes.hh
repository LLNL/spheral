//---------------------------------Spheral++----------------------------------//
// MortonOrderRedistributeNodes
//
// Attempt to redistribute nodes such that they are laid out in memory
// in a Morton ordering.  Note that this involves renumbering the nodes of 
// each NodeList, not just redistributing them between processors.
//
// Warren & Salmon (1995), Computer Physics Communications, 87, 266-290.
//
// Created by JMO, Tue Mar 25 14:19:18 PDT 2008
//----------------------------------------------------------------------------//
#ifndef MortonOrderRedistributeNodes_HH
#define MortonOrderRedistributeNodes_HH

#include "SpaceFillingCurveRedistributeNodes.hh"

namespace Spheral {
  template<typename Dimension> class DataBase;
}

namespace Spheral {

template<typename Dimension>
class MortonOrderRedistributeNodes: public SpaceFillingCurveRedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef KeyTraits::Key Key;

  // Constructors
  MortonOrderRedistributeNodes(const double dummy,
                               const double minNodesPerDomainFraction = 0.5,
                               const double maxNodesPerDomainFraction = 1.5,
                               const bool workBalance = true,
                               const bool localReorderOnly = false);

  // Destructor
  virtual ~MortonOrderRedistributeNodes();

  // Hash the positions.
  virtual
  FieldList<Dimension, Key> 
  computeHashedIndices(const DataBase<Dimension>& dataBase) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  // No copy or assignment operations.
  MortonOrderRedistributeNodes(const MortonOrderRedistributeNodes& nodes);
  MortonOrderRedistributeNodes& operator=(const MortonOrderRedistributeNodes& rhs);

};

}

#endif
