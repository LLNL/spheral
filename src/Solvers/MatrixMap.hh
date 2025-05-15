//---------------------------------Spheral++----------------------------------//
// MatrixMap
//
// Contains the information needed to initialize a matrix and linear solver.
// Assumes that the global indexing is sequential on one process and
// contiguous between processors.
//----------------------------------------------------------------------------//
#ifndef __Spheral_MatrixMap_hh__
#define __Spheral_MatrixMap_hh__

#include <vector>
#include "Utilities/DBC.hh"

namespace Spheral {

class MatrixMap {

public:
  // Constructor
  MatrixMap();
  
  // Destructor
  virtual ~MatrixMap() {}

  // Lower and upper global indices
  virtual int firstGlobalIndex() const = 0;
  virtual int lastGlobalIndex() const = 0;
  
  // Local and global number of elements
  virtual int numLocalElements() const = 0;
  virtual int numGlobalElements() const = 0;

  // Number of nonzero columns for each row
  virtual int numElementsPerRow(int localRowIndex) const = 0;
  virtual std::vector<int> numElementsPerRow() const;
  
  // Convert local to global index
  virtual int getGlobalIndex(int localIndex) const;

  // Get the nonzero column indices for a given row
  virtual void getColumnIndices(int localRowIndex,
                                std::vector<int> &globalColumnIndices) const {
    VERIFY2(false, "not implemented for this map type");
  }
  
  // Check class invariants
  virtual void checkClassInvariants() const = 0;
  
}; // end class MatrixMap
} // end namespace Spheral

#endif
