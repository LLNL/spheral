//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes -- (Re)domain decompose the nodes by using 
// a sort and divide algorithm.
//
// This is a template base class -- the actual dimension dependent objects
// are: 
//   SortAndDivideRedistributeNodes1d
//   SortAndDivideRedistributeNodes2d
//   SortAndDivideRedistributeNodes3d
//
// Created by JMO, Thu Dec 2 10:44:07 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_SortAndDivideRedistributeNodes_hh
#define Spheral_SortAndDivideRedistributeNodes_hh

#ifndef __GCCXML__
#include <list>
#include <vector>
#else
#include "fakestl.hh"
#endif

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "RedistributeNodes.hh"

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

template<typename Dimension>
class SortAndDivideRedistributeNodes: public RedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors
  SortAndDivideRedistributeNodes(const double HExtent);

  // Destructor
  virtual ~SortAndDivideRedistributeNodes();

  // Sort the given set of DomainNodes by a position component (where positionIndex
  // maps as
  //   0 ==> x,
  //   1 ==> y,
  //   2 ==> z),
  // and return the result as a list<DomainNode>.
  void sortByPositions(std::list<DomainNode<Dimension> >& domainNodes,
                       const int positionIndex) const;

  // Given a sorted list of domain nodes, pop nodes off of the front building
  // a return list until the requested work is met.
  std::list<DomainNode<Dimension> >
  popFrontNodes(std::list<DomainNode<Dimension> >& sortedCandidateNodes,
                const double targetDomainWork,
                const int positionIndex) const;

  // Compute the appropriate shape tensor for a set of domain nodes.
  typename SymTensor::EigenStructType shapeTensor(const std::vector<DomainNode<Dimension> >& domainNodes) const;

  // Apply the necessary rotation to the positions of the domain nodes to transform into
  // the primary frame of the given shape tensor.
  void rotateIntoShapeTensorFrame(const typename SymTensor::EigenStructType& shapeTensor,
                                  std::vector<DomainNode<Dimension> >& domainNodes) const;

  // Reduce a vector<DomainNode> to the given processor.
  std::vector<DomainNode<Dimension> > reduceDomainNodes(const std::vector<DomainNode<Dimension> >& nodes,
                                                        const int targetProc) const;

  // Broadcast a vector<DomainNode> from the given processor.
  std::vector<DomainNode<Dimension> > broadcastDomainNodes(const std::vector<DomainNode<Dimension> >& nodes,
                                                           const int sendProc) const;

  // Access the HExtent we're using.
  double Hextent() const;
  void Hextent(const double val);

protected:
  //--------------------------- Protected Interface ---------------------------//
#ifdef USE_MPI
  using RedistributeNodes<Dimension>::mCommunicator;
#endif

private:
  //--------------------------- Private Interface ---------------------------//
  // The cutoff radius in normalized space for nodes to interact.
  double mHextent;

  // No default constructor, copy, or assignment operations.
  SortAndDivideRedistributeNodes();
  SortAndDivideRedistributeNodes(const SortAndDivideRedistributeNodes&);
  SortAndDivideRedistributeNodes& operator=(const SortAndDivideRedistributeNodes&);

};

}
}

#ifndef __GCCXML__
#include "SortAndDivideRedistributeNodesInline.hh"
#endif

#else
// Forward declare the SortAndDivideRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    template<typename Dimension> class SortAndDivideRedistributeNodes;
  }
}

#endif
