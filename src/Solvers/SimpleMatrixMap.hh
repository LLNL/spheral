//---------------------------------Spheral++----------------------------------//
// SimpleMatrixMap
//
// Creates a map from the local rows and indices, including assigning
// global indices
//----------------------------------------------------------------------------//
#ifndef __Spheral_SimpleMatrixMap_hh__
#define __Spheral_SimpleMatrixMap_hh__

#include "MatrixMap.hh"

#include <memory>
#include <vector>

namespace Spheral {

class SimpleMatrixMap : public MatrixMap {
  
public:
  
  // Construtor.
  SimpleMatrixMap(int numLocalElements);

  // MatrixMap functions
  virtual int firstGlobalIndex() const override;
  virtual int lastGlobalIndex() const override;
  virtual int numLocalElements() const override;
  virtual int numGlobalElements() const override;
  virtual int numElementsPerRow(int localIndex) const override;
  // virtual std::vector<int> numElementsPerRow() const override;
  virtual int getGlobalIndex(int localIndex) const override;
  virtual void checkClassInvariants() const override;

  // Set row size
  virtual void setIndices(std::vector<std::vector<int>>& indices);

  // Get the nonzero column indices for a given row
  virtual void getColumnIndices(const int localRowIndex,
                                std::vector<int>& globalColumnIndices) const override;
  
  // Expose hidden base methods
  using MatrixMap::numElementsPerRow;

private:

  // Calculate global indices
  void initializeIndices();

  bool mIndicesSet;
  std::vector<std::vector<int>>* mIndices;
  
  int mFirstGlobalIndex;
  int mNumLocalElements;
  int mNumGlobalElements;
}; // end class MatrixData

} // end namespace Spheral

#endif
