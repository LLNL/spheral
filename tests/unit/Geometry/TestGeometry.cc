#include <iostream>
#include "GeomVector.hh"
#include "GeomTensor.hh"

void main () {

  // Test the GeomVector functions
  cerr << "############################## Vector ##############################" << endl;
  cerr << "GeomVector<3>(): " << GeomVector<3>() << endl;
  cerr << "GeomVector<3>(1.0): " << GeomVector<3>(1.0) << endl;

  GeomVector<3> r1(1.0, 2.0, 3.0);
  GeomVector<3> r2(2.0, 2.0, 2.0);
//   r1 = 1.0, 2.0, 3.0;
//   r2 = 2.0, 2.0, 2.0;
  cerr << "r1 = " << r1 << endl;
  cerr << "r2 = " << r2 << endl;

  GeomVector<3> r4(r1);
  cerr << "r4(r1) = " << r4 << endl;
  GeomVector<3> r3 = r2;
  cerr << "r3 = r2 = " << r3 << endl;
  double scalar = 5.0;
  r3 = scalar;
  cerr << "scalar = " << scalar << endl;
  cerr << "r3 = scalar = " << r3 << endl;
  cerr << "r1.x(), r1.y(), r1.z() = " << r1.x() << " "
       << r1.y() << " " << r1.z() << endl;
  r3.Zero();
  cerr << "r3.Zero(): " << r3 << endl;
  cerr << "-r1 = " << -r1 << endl;
  cerr << "r1 + r2 = " << r1 + r2 << endl;
  cerr << "r1 - r2 = " << r1 - r2 << endl;
  cerr << "r1*r2 = " << r1*r2 << endl;
  cerr << "r1 + scalar = " << r1 + scalar << endl;
  cerr << "r1 - scalar = " << r1 - scalar << endl;
  cerr << "r1*scalar = " << r1*scalar << endl;
  cerr << "r1/scalar = " << r1/scalar << endl;
  r3 = r1;
  r3 += r2;
  cerr << "r1 += r2: " << r3 << endl;
  r3 = r1;
  r3 -= r2;
  cerr << "r1 -= r2: " << r3 << endl;
  r3 = r1;
  r3 += scalar;
  cerr << "r1 += scalar: " << r3 << endl;
  r3 = r1;
  r3 -= scalar;
  cerr << "r1 -= scalar: " << r3 << endl;
  r3 = r1;
  r3 *= scalar;
  cerr << "r1 *= scalar: " << r3 << endl;
  r3 = r1;
  r3 /= scalar;
  cerr << "r1 /= scalar: " << r3 << endl;
  cerr << "scalar + r1 = " << scalar + r1 << endl;
  cerr << "scalar - r1 = " << scalar - r1 << endl;
  cerr << "scalar*r1 = " << scalar*r1 << endl;

  cerr << "r1 == r2: " << (r1 == r2) << endl;
  cerr << "r1 != r2: " << (r1 != r2) << endl;
  cerr << "r1 < r2: " << (r1 < r2) << endl;
  cerr << "r1 > r2: " << (r1 > r2) << endl;
  cerr << "r1 <= r2: " << (r1 <= r2) << endl;
  cerr << "r1 >= r2: " << (r1 >= r2) << endl;
  cerr << "r1 == r1: " << (r1 == r1) << endl;
  cerr << "r1 != r1: " << (r1 != r1) << endl;
  cerr << "r1 < r1: " << (r1 < r1) << endl;
  cerr << "r1 > r1: " << (r1 > r1) << endl;
  cerr << "r1 <= r1: " << (r1 <= r1) << endl;
  cerr << "r1 >= r1: " << (r1 >= r1) << endl;

  cerr << "r1.dot(r2) = " << r1.dot(r2) << endl;
  cerr << "r1.cross(r2) = " << r1.cross(r2) << endl;
  cerr << "r1.magnitude() = " << r1.magnitude() << endl;
  cerr << "r1.magnitude2() = " << r1.magnitude2() << endl;
  cerr << "r1.minElement() = " << r1.minElement() << endl;
  cerr << "r1.maxElement() = " << r1.maxElement() << endl;
  cerr << "r1.sumElements() = " << r1.sumElements() << endl;

  // Create a 3d full tensor.
  {
    cerr << "############################ Full Tensor ############################"
         << endl;
    cerr << "GeomTensor<3, FullTensor>(): " << endl
         << GeomTensor<3, FullTensor>() << endl;
    cerr << "GeomTensor<3, FullTensor>(1.0): " << endl
         << GeomTensor<3, FullTensor>(1.0) << endl;

    GeomTensor<3, FullTensor> T1(2.0, 1.0, 1.0,
                                 2.0, 1.0, 0.0,
                                 2.0, 0.0, 1.0);
    cerr << "T1 = " << endl << T1 << endl;

    GeomTensor<3, FullTensor> T2(1.0, 0.0, 0.0,
                                 2.0, 2.0, 0.0,
                                 3.0, 2.0, 3.0);
    cerr << "T2 = " << endl << T2 << endl;

    GeomTensor<3, FullTensor> T3;
    T3 = T1;
    cerr << "T3 = T1, T3 = " << endl << T3 << endl;
    GeomTensor<3, FullTensor> T4(T1);
    cerr << "T4(T1), T4 = " << endl << T4 << endl;
    T3 = scalar;
    cerr << "T3 = scalar, T3 = " << endl << T3 << endl;

    cerr << "Accessing T1 individual elements: "
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
    cerr << "T1.zz(5.0), T1 = " << endl << T3 << endl;

    cerr << "T1.getRow(1) = " << T1.getRow(1) << endl;
    cerr << "T1.getColumn(2) = " << T1.getColumn(2) << endl;

    T3 = T1;
    T3.setRow(1, GeomVector<3>(-1.0, -2.0, -3.0));
    cerr << "T1.setRow(1, (-1, -2, -3)), T1 = " << endl << T3 << endl;
    T3 = T1;
    T3.setColumn(2, GeomVector<3>(-5.0, -4.0, -3.0));
    cerr << "T1.setColumn(2, (-5, -4, -3)), T1 = " << endl << T3 << endl;

    T3 = T1;
    T3.Zero();
    cerr << "T1.Zero(), T1 = " << endl << T3 << endl;

    cerr << "-T1 = " << endl << -T1 << endl;
    cerr << "T1 + T2 = " << endl << T1 + T2 << endl;
    cerr << "T1 - T2 = " << endl << T1 - T2 << endl;
    cerr << "T1*T2 = " << endl << T1*T2 << endl;
    cerr << "T1*r1 = " << T1*r1 << endl;
    cerr << "T1 + scalar = " << endl << T1 + scalar << endl;
    cerr << "T1 - scalar = " << endl << T1 - scalar << endl;
    cerr << "T1*scalar = " << endl << T1*scalar << endl;
    cerr << "T1/scalar = " << endl << T1/scalar << endl;

    T3 = T1;
    T3 += T2;
    cerr << "T1 += T2, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 -= T2;
    cerr << "T1 -= T2, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 *= T2;
    cerr << "T1 *= T2, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 += scalar;
    cerr << "T1 += scalar, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 -= scalar;
    cerr << "T1 -= scalar, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 *= scalar;
    cerr << "T1 *= scalar, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 /= scalar;
    cerr << "T1 /= scalar, T1 = " << endl << T3 << endl;

    cerr << "T1 == T2: " << (T1 == T2) << endl;
    cerr << "T1 == T1: " << (T1 == T1) << endl;
    cerr << "T1 != T2: " << (T1 != T2) << endl;
    cerr << "T1 != T1: " << (T1 != T1) << endl;

    cerr << "T1.Symmetric() = " << endl << T1.Symmetric() << endl;
    cerr << "T1.SkewSymmetric() = " << endl << T1.SkewSymmetric() << endl;
    cerr << "T1.Transpose() = " << endl << T1.Transpose() << endl;
    cerr << "T1.Trace() = " << T1.Trace() << endl;
    cerr << "T1.Determinant() = " << T1.Determinant() << endl;
    cerr << "T1.dot(T2) = " << endl << T1.dot(T2) << endl;
    cerr << "T1.dot(r1) = " << T1.dot(r1) << endl;
    cerr << "T1.doubledot(T2) = " << T1.doubledot(T2) << endl;
  }

  // Create a 3d symmetric tensor.
  {
    cerr << "############################ Symmetric Tensor ############################"
         << endl;
    cerr << "GeomTensor<3, SymmetricTensor>(): " << endl
         << GeomTensor<3, SymmetricTensor>() << endl;
    cerr << "GeomTensor<3, SymmetricTensor>(1.0): " << endl
         << GeomTensor<3, SymmetricTensor>(1.0) << endl;

    GeomTensor<3, SymmetricTensor> T1(2.0, 0.0, 3.0,
                                      0.0, 1.0, 0.0,
                                      3.0, 0.0, 1.0);
    cerr << "T1 = " << endl << T1 << endl;

    GeomTensor<3, SymmetricTensor> T2(1.0, 3.0, 5.0,
                                      3.0, 2.0, 0.0,
                                      5.0, 0.0, 3.0);
    cerr << "T2 = " << endl << T2 << endl;

    GeomTensor<3, SymmetricTensor> T3;
    T3 = T1;
    cerr << "T3 = T1, T3 = " << endl << T3 << endl;
    GeomTensor<3, SymmetricTensor> T4(T1);
    cerr << "T4(T1), T4 = " << endl << T4 << endl;
    T3 = scalar;
    cerr << "T3 = scalar, T3 = " << endl << T3 << endl;

    cerr << "Accessing T1 individual elements: "
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
    cerr << "T1.zz(5.0), T1 = " << endl << T3 << endl;

    cerr << "T1.getRow(1) = " << T1.getRow(1) << endl;
    cerr << "T1.getColumn(2) = " << T1.getColumn(2) << endl;

    T3 = T1;
    T3.setRow(1, GeomVector<3>(-1.0, -2.0, -3.0));
    cerr << "T1.setRow(1, (-1, -2, -3)), T1 = " << endl << T3 << endl;
    T3 = T1;
    T3.setColumn(2, GeomVector<3>(-5.0, -4.0, -3.0));
    cerr << "T1.setColumn(2, (-5, -4, -3)), T1 = " << endl << T3 << endl;

    T3 = T1;
    T3.Zero();
    cerr << "T1.Zero(), T1 = " << endl << T3 << endl;

    cerr << "-T1 = " << endl << -T1 << endl;
    cerr << "T1 + T2 = " << endl << T1 + T2 << endl;
    cerr << "T1 - T2 = " << endl << T1 - T2 << endl;
    cerr << "T1*T2 = " << endl << T1*T2 << endl;
    cerr << "T1*r1 = " << T1*r1 << endl;
    cerr << "T1 + scalar = " << endl << T1 + scalar << endl;
    cerr << "T1 - scalar = " << endl << T1 - scalar << endl;
    cerr << "T1*scalar = " << endl << T1*scalar << endl;
    cerr << "T1/scalar = " << endl << T1/scalar << endl;

    T3 = T1;
    T3 += T2;
    cerr << "T1 += T2, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 -= T2;
    cerr << "T1 -= T2, T1 = " << endl << T3 << endl;
//     T3 = T1;
//     T3 *= T2;
//     cerr << "T1 *= T2, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 += scalar;
    cerr << "T1 += scalar, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 -= scalar;
    cerr << "T1 -= scalar, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 *= scalar;
    cerr << "T1 *= scalar, T1 = " << endl << T3 << endl;
    T3 = T1;
    T3 /= scalar;
    cerr << "T1 /= scalar, T1 = " << endl << T3 << endl;

    cerr << "T1 == T2: " << (T1 == T2) << endl;
    cerr << "T1 == T1: " << (T1 == T1) << endl;
    cerr << "T1 != T2: " << (T1 != T2) << endl;
    cerr << "T1 != T1: " << (T1 != T1) << endl;

    cerr << "T1.Symmetric() = " << endl << T1.Symmetric() << endl;
    cerr << "T1.SkewSymmetric() = " << endl << T1.SkewSymmetric() << endl;
    cerr << "T1.Transpose() = " << endl << T1.Transpose() << endl;
    cerr << "T1.Trace() = " << T1.Trace() << endl;
    cerr << "T1.Determinant() = " << T1.Determinant() << endl;
    cerr << "T1.dot(T2) = " << endl << T1.dot(T2) << endl;
    cerr << "T1.dot(r1) = " << T1.dot(r1) << endl;
    cerr << "T1.doubledot(T2) = " << T1.doubledot(T2) << endl;

    GeomTensor<3, FullTensor> T10(1, 2, 3,
                                  4, 5, 6,
                                  7, 8, 9);
    cerr << "T10 = " << endl << T10 << endl;
    cerr << "T1*T10 = " << endl << T1*T10 << endl;
    cerr << "T10*T1 = " << endl << T10*T1 << endl;

  }

  cerr << "Done." << endl;
}
  
