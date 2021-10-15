//------------------------------------------------------------------------------
// Standalone C++ executable to experiment with fuzzy comparisons.
//------------------------------------------------------------------------------

#include <iostream>
#include "../SpheralFunctions.hh"

using namespace Spheral;

int main() {

  const double xmin = numeric_limits<double>::min(),
               xmax = numeric_limits<double>::max(),
               eps = numeric_limits<double>::epsilon(),
               tol = 1.0e-10;
  auto x = xmin;
  while (x < xmax) {
    auto x2 = x + eps;
    auto x3 = 2.0*x;
    cout << " x <-> x + eps " << x << " : " << (x == x2) << " " << fuzzyEqual(x, x2, tol) << endl
         << " x <-> 2*x     " << x << " : " << (x == x3) << " " << fuzzyEqual(x, x3, tol) << endl << endl;
    x *= 10.0;
  }

  return 0;
}
