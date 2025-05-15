//---------------------------------Spheral++----------------------------------//
// SimpleMatrixData
//
// Creates a parallel matrix from the local rows and indices
//----------------------------------------------------------------------------//
#ifndef __Spheral_SimpleMatrixData_hh__
#define __Spheral_SimpleMatrixData_hh__

#include "MatrixData.hh"

#include <memory>
#include <vector>

namespace Spheral {

class SimpleMatrixMap;

class SimpleMatrixData : public MatrixData {
  
public:

  // Construtor.
  SimpleMatrixData(std::shared_ptr<SimpleMatrixMap> map);
  
  // Set the data
  virtual void setRowValues(int localRowIndex,
                            std::vector<int> const& globalColumnIndices,
                            std::vector<double> const& globalColumnValues);

  // MatrixData functions
  virtual std::shared_ptr<MatrixMap> getMap() const override;
  virtual void getRowValues(int localRowIndex,
                            std::vector<int>& globalColumnIndices,
                            std::vector<double>& globalColumnValues) const override;
  virtual void checkClassInvariants() const override;
  
  // Description of the data
  virtual std::string description() const override { return "SimpleMatrixData"; }

private:

  std::shared_ptr<SimpleMatrixMap> mMap;
  std::vector<std::vector<int> > mIndices;
  std::vector<std::vector<double> > mValues;
  
}; // end class MatrixData

} // end namespace Spheral

#endif
