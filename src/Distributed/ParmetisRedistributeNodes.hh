//---------------------------------Spheral++----------------------------------//
// ParmetisRedistributeNodes -- Redistribute nodes by calling the Parmetis
// library to determine the new optimal domain decomposition.
// http://www-users.cs.umn.edu/~karypis/metis/parmetis/
//
// Created by JMO, Wed Feb 12 17:29:09 PST 2003
//----------------------------------------------------------------------------//

#ifndef ParmetisRedistributeNodes_HH
#define ParmetisRedistributeNodes_HH

#ifndef __GCCXML__
#include <vector>
#include <map>
#else
#include "fakestl.hh"
#endif

#ifdef USE_MPI
#include "mpi.h"
extern "C" {
#include "parmetis.h"
}
#else
typedef int idx_t;
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
}

namespace Spheral {
namespace PartitionSpace {

template<typename Dimension>
class ParmetisRedistributeNodes: public RedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors
  ParmetisRedistributeNodes(double extent);

  // Destructor
  virtual ~ParmetisRedistributeNodes();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  // This is the method required of all descendent classes.
  virtual void redistributeNodes(DataBaseSpace::DataBase<Dimension>& dataBase,
                                 std::vector<BoundarySpace::Boundary<Dimension>*> boundaries =
                                 std::vector<BoundarySpace::Boundary<Dimension>*>());

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  // This is the method required of all descendent classes.
  void refineAndRedistributeNodes(DataBaseSpace::DataBase<Dimension>& dataBase,
                                  std::vector<BoundarySpace::Boundary<Dimension>*> boundaries =
                                  std::vector<BoundarySpace::Boundary<Dimension>*>());

  // The node extent we use when selecting neighbors.
  double normalizedNodeExtent() const;
  void setNormalizedNodeExtent(double extent);

private:
  //--------------------------- Private Interface ---------------------------//
  typedef std::map<int, std::vector<std::pair<int, double> > > ConnectivityType;

  // The cutoff radius in normalized space for nodes to interact.
  double mNormalizedNodeExtent;

  // Build the CSR adjacency information used by Parmetis.
  void buildCSRGraph(const DataBaseSpace::DataBase<Dimension>& dataBase,
                     const std::vector<DomainNode<Dimension> >& nodeDistribution,
                     const FieldSpace::FieldList<Dimension, int>& globalNodeIDs,
                     std::vector<idx_t>& xadj,
                     std::vector<idx_t>& adjacency,
                     std::vector<idx_t>& vtxdist,
                     std::vector<idx_t>& vweight,
                     std::vector<float>& xyz) const;

  // Cull the connectivity graph down to a dimension dependent number of connections per point.
  void cullConnectivity(ConnectivityType& neighbors,
                        const std::vector<DomainNode<Dimension> >& nodeDistribution) const;

  // Print the connectivity statistics.
  void printConnectivityStatistics(const ConnectivityType& neighbors) const;

  // Verify that the given connectivity data is minimally valid compared against a
  // domain decomposition.
  bool validConnectivity(const std::vector<DomainNode<Dimension> >& nodeDistribution,
                         const DataBaseSpace::DataBase<Dimension>& dataBase,
                         const ConnectivityType& neighbors,
                         const FieldSpace::FieldList<Dimension, int>& globalNodeIDs) const;

  // Verify that the given adjacency data is minimally valid compared against a
  // domain decomposition.
  bool validCSRGraph(const std::vector<DomainNode<Dimension> >& nodeDistribution,
                     const DataBaseSpace::DataBase<Dimension>& dataBase,
                     const std::vector<idx_t>& xadj,
                     const std::vector<idx_t>& adjacency,
                     const std::vector<idx_t>& vtxdist) const;

  // Helper method to check the given CSR graph lists a given node as connected
  // to another given node.
  bool verifyNeighborPresent(const int globalNeighborID,
                             const int globalID,
                             const std::vector<DomainNode<Dimension> >& nodeDistribution,
                             const std::vector<int>& xadj,
                             const std::vector<int>& adjacency,
                             const std::vector<int>& vtxdist) const;

  // Invert the given connectivity map.
  ConnectivityType inverseConnectivity(const ConnectivityType& neighbors) const;

  // Force the given connectivity to be symmetric. (ParMETIS requires it!)
  void enforceSymmetricConnectivity(ConnectivityType& neighbors) const;

  // Check the given connectivity is symmetric.
  bool symmetricConnectivity(const ConnectivityType& neighbors) const;

  // No copy or assignment operations.
  ParmetisRedistributeNodes(const ParmetisRedistributeNodes& nodes);
  ParmetisRedistributeNodes& operator=(const ParmetisRedistributeNodes& rhs);
};

}
}

#ifndef __GCCXML__
#include "ParmetisRedistributeNodesInline.hh"
#endif

#else
// Forward declare the ParmetisRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    template<typename Dimension> class ParmetisRedistributeNodes;
  }
}

#endif
