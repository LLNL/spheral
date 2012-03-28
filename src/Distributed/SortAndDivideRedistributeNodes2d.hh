//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes2d -- 2-D implementation of the sort and 
// divide algorithm for domain decomposition.
//
// Created by JMO, Sat Dec  4 16:28:49 PST 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_SortAndDivideRedistributeNodes2d_hh
#define Spheral_SortAndDivideRedistributeNodes2d_hh

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

class SortAndDivideRedistributeNodes2d: public SortAndDivideRedistributeNodes<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dim<2>::Scalar Scalar;
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::Tensor Tensor;
  typedef Dim<2>::SymTensor SymTensor;

  // Constructors
  SortAndDivideRedistributeNodes2d(const double HExtent);

  // Destructor
  virtual ~SortAndDivideRedistributeNodes2d();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  virtual void redistributeNodes(DataBaseSpace::DataBase<Dim<2> >& dataBase,
                                 std::vector<BoundarySpace::Boundary<Dim<2> >*> boundaries = std::vector<BoundarySpace::Boundary<Dim<2> >*>());

  // Compute the appropriate number of domains in each dimension for a given
  // shape tensor EigenStruct.
  std::vector<int> domainsPerChunk(const SymTensor::EigenStructType& shapeTensor) const;

private:
  //--------------------------- Private Interface ---------------------------//

#ifdef USE_MPI
  using RedistributeNodes<Dim<2> >::mCommunicator;
#endif

  // No default constructor, copy, or assignment operations.
  SortAndDivideRedistributeNodes2d();
  SortAndDivideRedistributeNodes2d(const SortAndDivideRedistributeNodes2d&);
  SortAndDivideRedistributeNodes2d& operator=(const SortAndDivideRedistributeNodes2d&);

};

}
}

#else
// Forward declare the SortAndDivideRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    class SortAndDivideRedistributeNodes2d;
  }
}

#endif
