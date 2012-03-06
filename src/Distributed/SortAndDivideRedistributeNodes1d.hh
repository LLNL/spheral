//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes1d -- 1-D implementation of the sort and 
// divide algorithm for domain decomposition.
//
// Created by JMO, Thu Dec 2 11:01:53 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_SortAndDivideRedistributeNodes1d_hh
#define Spheral_SortAndDivideRedistributeNodes1d_hh

#include "SortAndDivideRedistributeNodes.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace BoundarySpace {
    template<typename Dimension> class Boundary;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }
}

namespace Spheral {
namespace PartitionSpace {

class SortAndDivideRedistributeNodes1d: public SortAndDivideRedistributeNodes<Dim<1> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<1> Dimension;
  typedef Dim<1>::Scalar Scalar;
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::Tensor Tensor;
  typedef Dim<1>::SymTensor SymTensor;

  // Constructors
  SortAndDivideRedistributeNodes1d(const double HExtent);

  // Destructor
  virtual ~SortAndDivideRedistributeNodes1d();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  virtual void redistributeNodes(DataBaseSpace::DataBase<Dim<1> >& dataBase,
                                 std::vector<BoundarySpace::Boundary<Dim<1> >*> boundaries = std::vector<BoundarySpace::Boundary<Dim<1> >*>());


private:
  //--------------------------- Private Interface ---------------------------//

#ifdef USE_MPI
  using RedistributeNodes<Dim<1> >::mCommunicator;
#endif

  // No default constructor, copy, or assignment operations.
  SortAndDivideRedistributeNodes1d();
  SortAndDivideRedistributeNodes1d(const SortAndDivideRedistributeNodes1d&);
  SortAndDivideRedistributeNodes1d& operator=(const SortAndDivideRedistributeNodes1d&);

};

}
}

#else
// Forward declare the SortAndDivideRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    class SortAndDivideRedistributeNodes1d;
  }
}

#endif
