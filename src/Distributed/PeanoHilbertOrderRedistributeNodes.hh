//---------------------------------Spheral++----------------------------------//
// PeanoHilbertOrderRedistributeNodes
// Attempt to redistribute nodes such that they are laid out in memory
// in a PeanoHilbert ordering.  Note that this involves renumbering the nodes of 
// each NodeList, not just redistributing them between processors.
//
// Springel V. (2005), MNRAS
//
// Created by JMO, Tue Apr  8 14:25:18 PDT 2008
//----------------------------------------------------------------------------//
#ifndef PeanoHilbertOrderRedistributeNodes_HH
#define PeanoHilbertOrderRedistributeNodes_HH

#include "SpaceFillingCurveRedistributeNodes.hh"

namespace Spheral {
  template<typename Dimension> class DataBase;
}

namespace Spheral {

template<typename Dimension>
class PeanoHilbertOrderRedistributeNodes: public SpaceFillingCurveRedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename KeyTraits::Key Key;

  // Constructors
  PeanoHilbertOrderRedistributeNodes(const double dummy,
                                     const double minNodesPerDomainFraction = 0.5,
                                     const double maxNodesPerDomainFraction = 1.5,
                                     const bool workBalance = true,
                                     const bool localReorderOnly = false);

  // Destructor
  virtual ~PeanoHilbertOrderRedistributeNodes();

  // Hash the positions.
  virtual
  FieldList<Dimension, Key> 
  computeHashedIndices(const DataBase<Dimension>& dataBase) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  // No copy or assignment operations.
  PeanoHilbertOrderRedistributeNodes(const PeanoHilbertOrderRedistributeNodes& nodes);
  PeanoHilbertOrderRedistributeNodes& operator=(const PeanoHilbertOrderRedistributeNodes& rhs);

};

}

#endif
