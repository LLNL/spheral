//---------------------------------Spheral++----------------------------------//
// NestedGridRedistributeNodes -- (Re)domain decompose the nodes by using the
// information encoded in the NestedGridNeighbor algorithm.
//
// Created by JMO, Wed Nov 24 10:51:32 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_NestedGridRedistributeNodes_hh
#define Spheral_NestedGridRedistributeNodes_hh

#ifndef __GCCXML__
#include <set>
#include <vector>
#include <map>
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
  namespace NeighborSpace {
    template<typename Dimension> class GridCellIndex;
  }
}

namespace Spheral {
namespace PartitionSpace {

template<typename Dimension>
class NestedGridRedistributeNodes: public RedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef std::vector< std::map<NeighborSpace::GridCellIndex<Dimension>, int> > GridCellPopulationType;

  // Constructors
  NestedGridRedistributeNodes(const double HExtent);

  // Destructor
  virtual ~NestedGridRedistributeNodes();

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  // This is the method required of all RedistributeNodes classes.
  virtual void redistributeNodes(DataBaseSpace::DataBase<Dimension>& dataBase,
                                 std::vector<BoundarySpace::Boundary<Dimension>*> boundaries = std::vector<BoundarySpace::Boundary<Dimension>*>());

  // Find the global number of nodes on each grid level and in each grid cell.
  void computeGridCellPopulations(const DataBaseSpace::DataBase<Dimension>& dataBase,
                                  GridCellPopulationType& gridCellPopulations,
                                  std::vector<int>& gridLevelPopulations) const;

  // Find the grid level with the most nodes, and the minimum grid cell on that grid level.
  void findMaster(const GridCellPopulationType& gridCellPopulations,
                  const std::vector<int>& gridLevelPopulations,
                  NeighborSpace::GridCellIndex<Dimension>& gridCell,
                  int& gridLevel) const;

  // Count the number total number of nodes remaining in the given populations.
  int countRemainingNodes(const GridCellPopulationType& gridCellPopulations) const;
  int countRemainingNodes(const std::vector<int>& gridLevelPopulations) const;

  // Find the square "rind" of grid cells at the given distance from the given
  // grid cell on the given grid level.
  std::vector<NeighborSpace::GridCellIndex<Dimension> > gridCellRind(const NeighborSpace::GridCellIndex<Dimension>& gridCell,
                                                                     const int gridLevel,
                                                                     const int gridCellStep,
                                                                     const GridCellPopulationType& gridCellPopulations) const;

  // Helper method to compute the full set of GridCell's at a given step.
  std::set<NeighborSpace::GridCellIndex<Dimension> > computeGridCellRind(const NeighborSpace::GridCellIndex<Dimension>& gridCell,
                                                                         const int gridLevel) const;

  // Set the master grid cell info for all NodeLists in the DataBase.
  void setMasterNodeLists(DataBaseSpace::DataBase<Dimension>& dataBase,
                          const NeighborSpace::GridCellIndex<Dimension>& gridCell,
                          const int gridLevel) const;

  // Gather up the unassigned coarse neighbor nodes, filling in the global node indicies and work.
  void gatherAvailableCoarseNodes(const DataBaseSpace::DataBase<Dimension>& dataBase,
                                  const std::vector<DomainNode<Dimension> >& nodeDistribution,
                                  const FieldSpace::FieldList<Dimension, Scalar>& work,
                                  std::vector<int>& globalNodeIndicies,
                                  std::vector<Scalar>& globalNodeWork) const;

  // Assign from the given set of nodes to the given grid cell, until either the
  // desired work per domain is reached or we exhast the given set of nodes.
  bool assignNodesToDomain(const DataBaseSpace::DataBase<Dimension>& dataBase,
                           const std::vector<int>& globalNodeIndicies,
                           const std::vector<Scalar>& globalNodeWork,
                           const int currentDomainID,
                           const double targetWork,
                           double& curentWorkSum,
                           std::vector<DomainNode<Dimension> >& nodeDistribution,
                           GridCellPopulationType& gridCellPopulations,
                           std::vector<int>& gridLevelPopulations) const;

  // Access the HExtent we're using.
  double Hextent() const;
  void Hextent(const double val);

private:
  //--------------------------- Private Interface ---------------------------//
  // The cutoff radius in normalized space for nodes to interact.
  double mHextent;

  // No default constructor, copy, or assignment operations.
  NestedGridRedistributeNodes();
  NestedGridRedistributeNodes(const NestedGridRedistributeNodes&);
  NestedGridRedistributeNodes& operator=(const NestedGridRedistributeNodes&);

};

}
}

#ifndef __GCCXML__
#include "NestedGridRedistributeNodesInline.hh"
#endif

#else
// Forward declare the NestedGridRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    template<typename Dimension> class NestedGridRedistributeNodes;
  }
}

#endif
