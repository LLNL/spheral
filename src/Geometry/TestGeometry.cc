#include <iostream>
#include "GeomVector.hh"
#include "GeomTensor.hh"

int main () {

  // Test the GeomVector functions
  cout << "############################## Vector ##############################" << endl;
  cout << "GeomVector<3>(): " << GeomVector<3>() << endl;
  cout << "GeomVector<3>(1.0): " << GeomVector<3>(1.0) << endl;

  GeomVector<3> r1(1.0, 2.0, 3.0);
  GeomVector<3> r2(2.0, 2.0, 2.0);
//   r1 = 1.0, 2.0, 3.0;
//   r2 = 2.0, 2.0, 2.0;
  cout << "r1 = " << r1 << endl;
  cout << "r2 = " << r2 << endl;

  GeomVector<3> r3(r1);
  cout << "r3(r1) = " << r3 << endl;
  r3 = r2;
  cout << "r3 = r2 = " << r3 << endl;
  double scalar = 5.0;
  r3 = scalar;
  cout << "scalar = " << scalar << endl;
  cout << "r3 = scalar = " << r3 << endl;
  cout << "r1.x(), r1.y(), r1.z() = " << r1.x() << " "
       << r1.y() << " " << r1.z() << endl;
  r3.Zero();
  cout << "r3.Zero(): " << r3 << endl;
  cout << "-r1 = " << -r1 << endl;
  cout << "r1 + r2 = " << r1 + r2 << endl;
  cout << "r1 - r2 = " << r1 - r2 << endl;
  cout << "r1*r2 = " << r1*r2 << endl;
  cout << "r1 + scalar = " << r1 + scalar << endl;
  cout << "r1 - scalar = " << r1 - scalar << endl;
  cout << "r1*scalar = " << r1*scalar << endl;
  cout << "r1/scalar = " << r1/scalar << endl;
  r3 = r1;
  r3 += r2;
  cout << "r1 += r2: " << r3 << endl;
  r3 = r1;
  r3 -= r2;
  cout << "r1 -= r2: " << r3 << endl;
  r3 = r1;
  r3 += scalar;
  cout << "r1 += scalar: " << r3 << endl;
  r3 = r1;
  r3 -= scalar;
  cout << "r1 -= scalar: " << r3 << endl;
  r3 = r1;
  r3 *= scalar;
  cout << "r1 *= scalar: " << r3 << endl;
  r3 = r1;
  r3 /= scalar;
  cout << "r1 /= scalar: " << r3 << endl;
  cout << "scalar + r1 = " << scalar + r1 << endl;
  cout << "scalar - r1 = " << scalar - r1 << endl;
  cout << "scalar*r1 = " << scalar*r1 << endl;
  cout << "r1.dot(r2) = " << r1.dot(r2) << endl;
  cout << "r1.cross(r2) = " << r1.cross(r2) << endl;
  cout << "r1.magnitude() = " << r1.magnitude() << endl;
  cout << "r1.magnitude2() = " << r1.magnitude2() << endl;
  // Create a 3d tensor.
  cout << "############################## Tensor ##############################" << endl;
  cout << "GeomTensor<3>(): " << GeomTensor<3>() << endl;
  cout << "GeomTensor<3>(1.0): " << GeomTensor<3>(1.0) << endl;

  GeomTensor<3> T1(2.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
  cout << "T1 = " << T1 << endl;

  GeomTensor<3> T2(1.0, 0.0, 0.0,
                   0.0, 2.0, 0.0,
                   0.0, 0.0, 3.0);
  cout << "T2 = " << T2 << endl;

  GeomTensor<3> T3;
  T3 = T1;
  cout << "T3 = T1, T3 = " << T3 << endl;
  GeomTensor<3> T4(T1);
  cout << "T4(T1), T4 = " << T4 << endl;
  T3 = scalar;
  cout << "T3 = scalar, T3 = " << T3 << endl;

  cout << "Accessing T1 individual elements: "
       << T1.xx() << " "
       << T1.xy() << " "
       << T1.xz() << " "
       << T1.yx() << " "
       << T1.yy() << " "
       << T1.yz() << " "
       << T1.zx() << " "
       << T1.zy() << " "
       << T1.zz() << endl;

  T3 = T1;
  T3.zz(5.0);
  cout << "T1.zz(5.0), T1 = " << T3 << endl;

  cout << "T1.row(1) = " << T1.row(1) << endl;
  cout << "T1.column(2) = " << T1.column(2) << endl;

  T3 = T1;
  T3.row(1, GeomVector<3>(1.0, 2.0, 3.0));
  cout << "T1.row(1, (1, 2, 3)), T1 = " << T3 << endl;
  T3 = T1;
  T3.column(2, GeomVector<3>(5.0, 4.0, 3.0));
  cout << "T1.column(2, (5, 4, 3)), T1 = " << T3 << endl;

  T3 = T1;
  T3.Zero();
  cout << "T1.Zero(), T1 = " << T3 << endl;

  cout << "-T1 = " << -T1 << endl;
  cout << "T1 + T2 = " << T1 + T2 << endl;
  cout << "T1 - T2 = " << T1 - T2 << endl;
  cout << "T1*T2 = " << T1*T2 << endl;
  cout << "T1*r1 = " << T1*r1 << endl;
  cout << "T1 + scalar = " << T1 + scalar << endl;
  cout << "T1 - scalar = " << T1 - scalar << endl;
  cout << "T1*scalar = " << T1*scalar << endl;
  cout << "T1/scalar = " << T1/scalar << endl;

  T3 = T1;
  T3 += T2;
  cout << "T1 += T2, T1 = " << T3 << endl;
  T3 = T1;
  T3 -= T2;
  cout << "T1 -= T2, T1 = " << T3 << endl;
  T3 = T1;
  T3 *= T2;
  cout << "T1 *= T2, T1 = " << T3 << endl;
  T3 = T1;
  T3 += scalar;
  cout << "T1 += scalar, T1 = " << T3 << endl;
  T3 = T1;
  T3 -= scalar;
  cout << "T1 -= scalar, T1 = " << T3 << endl;
  T3 = T1;
  T3 *= scalar;
  cout << "T1 *= scalar, T1 = " << T3 << endl;
  T3 = T1;
  T3 /= scalar;
  cout << "T1 /= scalar, T1 = " << T3 << endl;

  cout << "T1.Symmetric() = " << T1.Symmetric() << endl;
  cout << "T1.SkewSymmetric() = " << T1.SkewSymmetric() << endl;
  cout << "T1.Transpose() = " << T1.Transpose() << endl;
  cout << "T1.Trace() = " << T1.Trace() << endl;
  cout << "T1.Determinant() = " << T1.Determinant() << endl;
  cout << "T1.dot(T2) = " << T1.dot(T2) << endl;
  cout << "T1.dot(r1) = " << T1.dot(r1) << endl;
  cout << "T1.doubledot(T2) = " << T1.doubledot(T2) << endl;

  return 0;
}
  
