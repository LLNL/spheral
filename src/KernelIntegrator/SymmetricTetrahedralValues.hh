//---------------------------------Spheral++----------------------------------//
// SymmetricTetrahedralValues
//
// Symmetric tetrahedreal quadrature with positive weights
// https://doi.org/10.1016/j.camwa.2015.03.017
//----------------------------------------------------------------------------//
#ifndef __Spheral_SymmetricTetrahedralValues_hh__
#define __Spheral_SymmetricTetrahedralValues_hh__

#include <vector>

#include "Geometry/Dimension.hh"

namespace Spheral {

struct SymmetricTetrahedralValues {
  static const std::vector<double> values1;  // order 1
  static const std::vector<double> values4;  // order 2
  static const std::vector<double> values8;  // order 3
  static const std::vector<double> values14; // order 5
  static const std::vector<double> values24; // order 6
  static const std::vector<double> values35; // order 7
  static const std::vector<double> values46; // order 8
  static const std::vector<double> values59; // order 9
  static const std::vector<double> values81; // order 10

  static int numOrdinatesForOrder(const int order);
  static int orderForNumOrdinates(const int numOrdinates);
  
  // Get values from above list
  static const std::vector<double>& getValues(const int numOrdinates);
  
  // Get quadrature for x+y\in(0,1)
  static void getQuadrature(const int numOrdinates,
                            std::vector<double>& weights,
                            std::vector<Dim<3>::Vector>& ordinates);
};

} // end namespace Spheral

#endif
