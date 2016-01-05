#include <iostream>

#include "mtl/mtl.h"
#include "mtl/utils.h"
#include "mtl/matrix.h"
#include "mtl/linalg_vec.h"
using namespace mtl;

int main () {
  double r1data[] = {1.0, 1.0, 1.0};
  dense1D<double> r1(3);
  mtl::set(r1, 1.0);
  dense1D<double> r2(r1);

  cout << "r1 = " << &r1 << " ";
  print_vector(r1);
  cout << "r2 = " << &r2 << " ";
  print_vector(r2);

  double scalar = 5.0;
  mtl::scale(r2, scalar);

  cout << "r1 = " << &r1 << " ";
  print_vector(r1);
  cout << "r2 = " << &r2 << " ";
  print_vector(r2);

  return 0;
}  
