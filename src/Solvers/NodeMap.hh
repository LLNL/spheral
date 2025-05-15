//---------------------------------Spheral++----------------------------------//
// SpheralNodeMap
//
// Creates a matrix map based on the connectivity of the Spheral nodes.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeMap_hh__
#define __Spheral_NodeMap_hh__

#include <utility>
#include <vector>
#include "Field/FieldList.hh"
#include "KernelIntegrator/FlatConnectivity.hh"
#include "Registrar/BoundaryDependent.hh"
#include "Registrar/ConnectivityDependent.hh"
#include "Registrar/UpdateDependent.hh"
#include "MatrixMap.hh"

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class CombinedRegistrar;

template<typename Dimension>
class NodeMap : public MatrixMap {
  
public:
  // Constructor
  NodeMap(const FlatConnectivity<Dimension>& flatConnectivity);

  // Destructor
  virtual ~NodeMap() {}

  // Lower and upper global indices
  virtual int firstGlobalIndex() const override;
  virtual int lastGlobalIndex() const override;
  
  // Local and global number of elements
  virtual int numLocalElements() const override;
  virtual int numGlobalElements() const override;

  // Number of nonzero columns for this row
  virtual int numElementsPerRow(int localRowIndex) const override;
  virtual int numElementsPerRow(int nodeListIndex,
                                int nodeIndex) const;

  // Check class invariants
  virtual void checkClassInvariants() const override;

  // Get the nonzero column indices for a given row
  virtual void getColumnIndices(const int localRowIndex,
                                std::vector<int>& globalColumnIndices) const override;
  
  // Get local index from NodeList and Node indices
  virtual int getLocalIndex(int nodeListIndex,
                            int nodeIndex) const;

  // Get NodeList and Node indices from local index
  virtual std::pair<int, int> getNodeIndex(int localIndex) const;

  // Get global index from NodeList and Node indices
  virtual int getGlobalIndex(int localIndex) const override;
  virtual int getGlobalIndex(int nodeListIndex,
                             int nodeIndex) const;
  
  // Get boundary information
  virtual bool isConstantBoundaryNode(int nodeListIndex,
                                      int nodeIndex) const;

  // Return connectivity
  virtual const FlatConnectivity<Dimension>& flatConnectivity() const { return mConnectivity; }
  
  // virtual const std::vector<int>& constantBoundaryNodes(int nodeListIndex) const;
  
  // Expose hidden base methods
  using MatrixMap::numElementsPerRow;

protected:

  // Data
  const FlatConnectivity<Dimension>& mConnectivity;
  
}; // end class NodeMap

} // end namespace Spheral

#endif
