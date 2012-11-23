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
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {
namespace PartitionSpace {

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
  FieldSpace::FieldList<Dimension, Key> 
  computeHashedIndicies(const DataBaseSpace::DataBase<Dimension>& dataBase) const;

protected:
  //--------------------------- Protected Interface ---------------------------//
#ifdef USE_MPI
  using RedistributeNodes<Dimension>::mCommunicator;
#endif

private:
  //--------------------------- Private Interface ---------------------------//
  // No copy or assignment operations.
  MortonOrderRedistributeNodes(const MortonOrderRedistributeNodes& nodes);
  MortonOrderRedistributeNodes& operator=(const MortonOrderRedistributeNodes& rhs);

};

}
}

#else
// Forward declare the MortonOrderRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    template<typename Dimension> class MortonOrderRedistributeNodes;
  }
}

#endif
