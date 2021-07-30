//---------------------------------Spheral++----------------------------------//
// NeighborNodeList -- A NodeList with neighbors and connectivity
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborNodeList__
#define __Spheral_NeighborNodeList__

#include "NodeList.hh"

#include <float.h>
#include <string>

namespace Spheral {

template<typename Dimension> class NodeIteratorBase;
template<typename Dimension> class AllNodeIterator;
template<typename Dimension> class InternalNodeIterator;
template<typename Dimension> class GhostNodeIterator;
template<typename Dimension> class MasterNodeIterator;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class Neighbor;
template<typename Dimension> class ConnectivityMap;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class EquationOfState;
template<typename Dimension> class TableKernel;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class DataBase;
class FileIO;

template<typename Dimension>
class NeighborNodeList: public NodeList<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Constructors.
  NeighborNodeList(std::string name,
                   const int numInternal,
                   const int numGhost,
                   const Scalar hmin = 1.0e-20,
                   const Scalar hmax = 1.0e20,
                   const Scalar hminratio = 0.1,
                   const Scalar nPerh = 2.01,
                   const int maxNumNeighbors = 500);

  // Destructor
  virtual ~NeighborNodeList();

  // Node iterators for neighbor nodes
  MasterNodeIterator<Dimension> masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const;
  MasterNodeIterator<Dimension> masterNodeEnd() const;
          
  CoarseNodeIterator<Dimension> coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const;
  CoarseNodeIterator<Dimension> coarseNodeEnd() const;

  RefineNodeIterator<Dimension> refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const;
  RefineNodeIterator<Dimension> refineNodeEnd() const;

  // The maximum number of neighbors we want to have (for calculating the ideal H).
  unsigned maxNumNeighbors() const;
  void maxNumNeighbors(unsigned val);

  // Access the neighbor object.
  Neighbor<Dimension>& neighbor() const;

  void registerNeighbor(Neighbor<Dimension>& neighbor);
  void unregisterNeighbor();

  //****************************************************************************
  // Methods required for restarting.
  // Dump and restore the NodeList state.
  virtual std::string label() const { return "NeighborNodeList"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
   // Neighbor information
  unsigned mMaxNumNeighbors;
  Neighbor<Dimension>* mNeighborPtr;
   
  // No default constructor or copying.
  NeighborNodeList();
  NeighborNodeList(const NeighborNodeList& nodes);
  NeighborNodeList& operator=(const NeighborNodeList& rhs);
};

}

#include "NeighborNodeListInline.hh"

#else

// Forward declaration
namespace Spheral {
  template<typename Dimension> class NeighborNodeList;
}

#endif

