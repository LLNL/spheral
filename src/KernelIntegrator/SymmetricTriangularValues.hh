//---------------------------------Spheral++----------------------------------//
// SymmetricTriangularValues
//
// Symmetric triangular quadrature with positive weights
// https://doi.org/10.1016/j.camwa.2015.03.017
//----------------------------------------------------------------------------//
#ifndef __Spheral_SymmetricTriangularValues_hh__
#define __Spheral_SymmetricTriangularValues_hh__

#include <vector>

#include "Geometry/Dimension.hh"

namespace Spheral {

struct SymmetricTriangularValues {
  static const std::vector<double> values1;  // order 1
  static const std::vector<double> values3;  // order 2
  static const std::vector<double> values6;  // order 4
  static const std::vector<double> values7;  // order 5
  static const std::vector<double> values12; // order 6
  static const std::vector<double> values15; // order 7
  static const std::vector<double> values16; // order 8 
  static const std::vector<double> values19; // order 9
  static const std::vector<double> values25; // order 10

  // Get order and number of ordinates in terms of one another
  static int numOrdinatesForOrder(const int order);
  static int orderForNumOrdinates(const int numOrdinates);
  
  // Get values from above list
  static const std::vector<double>& getValues(const int numOrdinates);
  
  // Get quadrature for x+y\in(0,1)
  static void getQuadrature(const int numOrdinates,
                            std::vector<double>& weights,
                            std::vector<Dim<2>::Vector>& ordinates);
  static void getQuadrature(const int numOrdinates,
                            std::vector<double>& weights,
                            std::vector<Dim<3>::Vector>& ordinates);
};

} // end namespace Spheral

#endif
