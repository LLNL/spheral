//---------------------------------Spheral++----------------------------------//
// DistributeByXPosition -- Redistribute nodes by sorting their positions
// in x coordinate.  Really only useful in 1-D, as a test.
//
// Created by JMO, Mon Feb 10 16:12:47 PST 2003
//----------------------------------------------------------------------------//

#ifndef DistributeByXPosition_HH
#define DistributeByXPosition_HH

#include "RedistributeNodes.hh"

namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
}

namespace Spheral {
namespace PartitionSpace {

template<typename Dimension>
class DistributeByXPosition: public RedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors
  DistributeByXPosition();

  // Destructor
  virtual ~DistributeByXPosition();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  // This is the method required of all descendent classes.
  virtual void redistributeNodes(DataBaseSpace::DataBase<Dimension>& dataBase,
                                 std::vector<BoundarySpace::Boundary<Dimension>*> boundaries = std::vector<BoundarySpace::Boundary<Dimension>*>());

private:
  //--------------------------- Private Interface ---------------------------//
  // No copy or assignment operations.
  DistributeByXPosition(const DistributeByXPosition& nodes);
  DistributeByXPosition& operator=(const DistributeByXPosition& rhs);

#ifdef USE_MPI
  using RedistributeNodes<Dimension>::mCommunicator;
#endif

};

}
}

#else
// Forward declare the DistributeByXPosition class.
namespace Spheral {
  namespace PartitionSpace {
    template<typename Dimension> class DistributeByXPosition;
  }
}

#endif
