//---------------------------------Spheral++----------------------------------//
// MatrixData
//
// Represents the data needed to define an explicit matrix for multiplication
// or inversion.
//----------------------------------------------------------------------------//
#ifndef __Spheral_MatrixData_hh__
#define __Spheral_MatrixData_hh__

#include <memory>
#include <vector>
#include <string>

namespace Spheral {

class MatrixMap;

class MatrixData : public std::enable_shared_from_this<MatrixData> {
  
public:
  // Construtor.
  MatrixData();
  
  // Get the matrix map
  virtual std::shared_ptr<MatrixMap> getMap() const = 0;

  // Return the nonzero indices and values of the matrix.
  virtual void getRowValues(int localRowIndex,
                            std::vector<int>& globalColumnIndices,
                            std::vector<double>& globalColumnValues) const = 0;

  // Optionally allows for setting of initial guess to something other than the RHS
  virtual bool initialGuessAvailable() const;
  virtual void getInitialGuess(std::vector<double>& data) const;
  
  // Check class invariants.
  virtual void checkClassInvariants() const = 0;

  // Description of the data
  virtual std::string description() const { return "MatrixData"; }
  
protected:
  // Get shared_ptr for this object
  std::shared_ptr<MatrixData> getPointer();
  
}; // end class MatrixData

} // end namespace Spheral

#include "MatrixDataInline.hh"

#endif
