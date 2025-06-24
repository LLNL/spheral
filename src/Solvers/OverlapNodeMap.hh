//---------------------------------Spheral++----------------------------------//
// OverlapNodeMap
//
// Node map for overlap connectivity
//----------------------------------------------------------------------------//
#ifndef __Spheral_OverlapNodeMap_hh__
#define __Spheral_OverlapNodeMap_hh__

#include <vector>
#include "MatrixMap.hh"
#include "NodeMap.hh"
#include "KernelIntegrator/FlatConnectivity.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

template<typename Dimension>
class OverlapNodeMap : public NodeMap<Dimension> {

public:
  // Constructor
  OverlapNodeMap(const FlatConnectivity<Dimension>& connectivity);


  // Number of nonzero columns for this row
  virtual int numElementsPerRow(int localRowIndex) const override;
  virtual int numElementsPerRow(int nodeListIndex,
                                int nodeIndex) const override;
  
  // Get the nonzero column indices for a given row
  virtual void getColumnIndices(int localRowIndex,
                                std::vector<int> &globalColumnIndices) const override;

  // Expose hidden methods
  using MatrixMap::numElementsPerRow;
  
}; // end class OverlapNodeMap

} // end namespace Spheral

#endif
