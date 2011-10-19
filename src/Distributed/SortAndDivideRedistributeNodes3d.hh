//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes3d -- 3-D implementation of the sort and 
// divide algorithm for domain decomposition.
//
// Created by JMO, Sun Dec  5 21:44:28 PST 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_SortAndDivideRedistributeNodes3d_hh
#define Spheral_SortAndDivideRedistributeNodes3d_hh

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

class SortAndDivideRedistributeNodes3d: public SortAndDivideRedistributeNodes<Dim<3> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3> Dimension;
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;

  // Constructors
  SortAndDivideRedistributeNodes3d(const double HExtent);

  // Destructor
  virtual ~SortAndDivideRedistributeNodes3d();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  virtual void redistributeNodes(DataBaseSpace::DataBase<Dim<3> >& dataBase,
                                 std::vector<BoundarySpace::Boundary<Dim<3> >*> boundaries = std::vector<BoundarySpace::Boundary<Dim<3> >*>());

  // Compute the appropriate number of domains in each dimension for a given
  // shape tensor EigenStruct.
  std::vector< std::vector<int> > domainsPerChunk(const SymTensor::EigenStructType& shapeTensor) const;

private:
  //--------------------------- Private Interface ---------------------------//

#ifdef USE_MPI
  using RedistributeNodes<Dim<3> >::mCommunicator;
#endif

  // No default constructor, copy, or assignment operations.
  SortAndDivideRedistributeNodes3d();
  SortAndDivideRedistributeNodes3d(const SortAndDivideRedistributeNodes3d&);
  SortAndDivideRedistributeNodes3d& operator=(const SortAndDivideRedistributeNodes3d&);

};

}
}

#else
// Forward declare the SortAndDivideRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    class SortAndDivideRedistributeNodes3d;
  }
}

#endif
