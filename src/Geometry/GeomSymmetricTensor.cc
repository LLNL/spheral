//---------------------------------Spheral++----------------------------------//
// GeomSymmetricTensor -- the symmetric tensor class
//----------------------------------------------------------------------------//
#include <cmath>
#include <limits.h>
#include <float.h>
#include <vector>
using namespace std;

#include "GeomSymmetricTensor.hh"
#include "EigenStruct.hh"
#include "buildEigenVector.hh"
#include "findEigenValues3.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "Infrastructure/SpheralError.hh"
#include "Utilities/rotationMatrix.hh"

#include "Jacobi2.hh"

#include <Eigen/Dense>

namespace Spheral {

using std::abs;

using Geometry::buildUniqueEigenVector;
using Geometry::jacobiDiagonalize;

//------------------------------------------------------------------------------
// Return the eigen values and eigen vectors of a symmetric tensor
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// 3-D.
template<>
EigenStruct<3>
GeomSymmetricTensor<3>::eigenVectors() const {

  // Some useful typedefs.
  typedef GeomVector<3> Vector;
  typedef GeomTensor<3> Tensor;
  typedef GeomSymmetricTensor<3> SymTensor;

  // Tolerances for fuzzy math.
  const double degenerate = 1.0e-20;
  const double tolerance = 5.0e-5;

  // Prepare the result.
  EigenStruct<3> result;

  // Create a scaled version of this tensor, with all elements in the range [-1,1].
  const double fscale = max(1.0, this->maxAbsElement());
  CHECK(fscale >= 1.0);
  const double fscalei = 1.0/fscale;
  const SymTensor A = (*this)*fscalei;

// #ifdef USEJACOBI

//   // Use the Jacobi iterative diagonalization method to determine
//   // the eigen values/vectors.
//   const int nrot = jacobiDiagonalize<Dim<3> >(A,
//                                               result.eigenVectors,
//                                               result.eigenValues);
//   result.eigenValues *= fscale;

// #elif USEEIGEN

  // Use the Eigen library to determine the eigen values/vectors.
  {
    Eigen::Matrix3d B;
    B << 
      A.xx(), A.xy(), A.xz(),
      A.yx(), A.yy(), A.yz(),
      A.zx(), A.zy(), A.zz();
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(B);
    const Eigen::Vector3d& Bvals = eigensolver.eigenvalues();
    const Eigen::Matrix3d& Bvecs = eigensolver.eigenvectors();
    result.eigenValues = Vector(Bvals(0), Bvals(1), Bvals(2)) * fscale;
    result.eigenVectors = Tensor(Bvecs(0,0), Bvecs(0,1), Bvecs(0,2),
                                 Bvecs(1,0), Bvecs(1,1), Bvecs(1,2),
                                 Bvecs(2,0), Bvecs(2,1), Bvecs(2,2));
  }

// #else

//   // Compute the scaled eigen-values, and sort them.
//   Vector lambdaVec = A.eigenValues();
//   sort(lambdaVec.begin(), lambdaVec.end());
//   CHECK(lambdaVec.x() <= lambdaVec.y() and
//         lambdaVec.y() <= lambdaVec.z());

//   // Assign the true eigen-values in the result.
//   result.eigenValues = fscale*lambdaVec;
//   result.eigenVectors = SymTensor::one;

//   // If any of the eigen-values result in a tensor that is not positive-rank 
//   // (all zero elements), we assume the eigen-values are equal and punt
//   // with the identity tensor for the eigen-vectors.
//   // We simultaneously compute the row containing the maximum absolute value 
//   // element for each eigen-value.
//   bool punt = false;
//   double maxEVelement = -1.0;
//   Vector maxEVrow;
//   int iFirst = -1;
//   for (int ivalue = 0; ivalue != 3; ++ivalue) {
//     const SymTensor M = A - lambdaVec(ivalue)*SymTensor::one;
//     if (M.maxAbsElement() < degenerate) punt = true;
//     for (int irow = 0; irow != 3; ++irow) {
//       const Vector Mvec = M.getRow(irow);
//       const double thpt = Mvec.maxAbsElement();
//       if (thpt > maxEVelement) {
//         maxEVelement = thpt;
//         maxEVrow = Mvec;
//         iFirst = ivalue;
//       }
//     }
//   }

//   // If we found an all zero M (= A - lambda*I) matrix, we punt and accept the identity
//   // tensor as our eigen-vectors.  Otherwise, continue the compuation.
//   if (!punt) {
//     CHECK(iFirst >= 0 and iFirst < 3);

//     // Select the ordering we'll go through the eigen-values in, starting
//     // with the row with the largest absolute value element.
//     const int iSecond = (iFirst + 1) % 3;
//     const int iThird = (iSecond + 1) % 3;
//     CHECK(iFirst + iSecond + iThird == 3);

//     // We need two orthogonal unit vectors in the plane perpendicular to
//     // the maximum row selected previously.  We can do this by finding the
//     // rotational transformation wherein x' axis is aligned with this row, and 
//     // taking our two vectors as the other two rows of this transform.
//     const Vector R = maxEVrow.unitVector();
//     const Tensor Tr = rotationMatrix(R);
//     const Vector U0 = Tr.getRow(1);
//     const Vector U1 = Tr.getRow(2);
    
//     // Now we can compute the eigen-vector corresponding the first eigen-value
//     // selected previously.
//     const Vector V0 = buildUniqueEigenVector(A, 
//                                              lambdaVec(iFirst),
//                                              U0,
//                                              U1);
//     result.eigenVectors.setColumn(iFirst, V0);

//     // Now we know the remaining eigen-vectors are in the plane perpendicular to
//     // V0.  We know R is in that plane, and so is R x V0.  With that knowledge
//     // we can basically repeat the same procedure for the next eigen-vector.
//     Vector S = R.cross(V0);
//     CHECK(fuzzyEqual(S.magnitude2(), 1.0, tolerance));
//     const Vector V1 = buildUniqueEigenVector(A,
//                                              lambdaVec(iSecond),
//                                              R,
//                                              S);
//     result.eigenVectors.setColumn(iSecond, V1);
    
//     // The last eigen-vector is orthogonal to the first two, so we can find it
//     // simply by taking the cross-product of the previous eigen-vectors.
//     const Vector V2 = V0.cross(V1);
//     CHECK(fuzzyEqual(V2.magnitude2(), 1.0, tolerance));
//     CHECK(fuzzyEqual(((A - lambdaVec(iThird)*SymTensor::one)*V2).maxAbsElement(), 0.0, tolerance));
//     result.eigenVectors.setColumn(iThird, V2);
//   }

// #endif

  BEGIN_CONTRACT_SCOPE;
  // Check the result.
  const double lambda1 = result.eigenValues.x();
  const double lambda2 = result.eigenValues.y();
  const double lambda3 = result.eigenValues.z();
  const Vector v1 = result.eigenVectors.getColumn(0);
  const Vector v2 = result.eigenVectors.getColumn(1);
  const Vector v3 = result.eigenVectors.getColumn(2);
  ENSURE(fuzzyEqual(v1.dot(v2), 0.0, tolerance) and 
         fuzzyEqual(v1.dot(v3), 0.0, tolerance) and 
         fuzzyEqual(v2.dot(v3), 0.0, tolerance));
  ENSURE(fuzzyEqual(v1.magnitude2(), 1.0, tolerance) and
         fuzzyEqual(v2.magnitude2(), 1.0, tolerance) and
         fuzzyEqual(v3.magnitude2(), 1.0, tolerance));
  const double tol = tolerance*max(1.0, this->maxAbsElement());
  ENSURE2(fuzzyEqual((SymTensor(xx() - lambda1, xy(), xz(),
                                yx(), yy() - lambda1, yz(),
                                zx(), zy(), zz() - lambda1)*v1).maxAbsElement(), 0.0, tol),
          *this << " " << lambda1 << " " << v1 << " " << tol << " "
          << SymTensor(xx() - lambda1, xy(), xz(),
                       yx(), yy() - lambda1, yz(),
                       zx(), zy(), zz() - lambda1)*v1);
  ENSURE(fuzzyEqual((SymTensor(xx() - lambda2, xy(), xz(),
                               yx(), yy() - lambda2, yz(),
                               zx(), zy(), zz() - lambda2)*v2).maxAbsElement(), 0.0, tol));
  ENSURE(fuzzyEqual((SymTensor(xx() - lambda3, xy(), xz(),
                               yx(), yy() - lambda3, yz(),
                               zx(), zy(), zz() - lambda3)*v3).maxAbsElement(), 0.0, tol));
  ENSURE(fuzzyEqual(abs(result.eigenVectors.Determinant()), 1.0, tolerance));
  END_CONTRACT_SCOPE;

  return result;
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class GeomSymmetricTensor<1>;
template class GeomSymmetricTensor<2>;
template class GeomSymmetricTensor<3>;

//------------------------------------------------------------------------------
// Set the static variables.
//------------------------------------------------------------------------------
template<> const int GeomSymmetricTensor<1>::nDimensions = 1;
template<> const GeomSymmetricTensor<1> GeomSymmetricTensor<1>::zero = GeomSymmetricTensor<1>(0.0);
template<> const GeomSymmetricTensor<1> GeomSymmetricTensor<1>::one = GeomSymmetricTensor<1>(1.0);
template<> const GeomSymmetricTensor<1>::size_type GeomSymmetricTensor<1>::mNumElements = 1;

template<> const int GeomSymmetricTensor<2>::nDimensions = 2;
template<> const GeomSymmetricTensor<2> GeomSymmetricTensor<2>::zero = GeomSymmetricTensor<2>(0.0, 0.0,
                                                                                              0.0, 0.0);
template<> const GeomSymmetricTensor<2> GeomSymmetricTensor<2>::one = GeomSymmetricTensor<2>(1.0, 0.0,
                                                                                             0.0, 1.0);
template<> const GeomSymmetricTensor<2>::size_type GeomSymmetricTensor<2>::mNumElements = 3;

template<> const int GeomSymmetricTensor<3>::nDimensions = 3;
template<> const GeomSymmetricTensor<3> GeomSymmetricTensor<3>::zero = GeomSymmetricTensor<3>(0.0, 0.0, 0.0,
                                                                                              0.0, 0.0, 0.0,
                                                                                              0.0, 0.0, 0.0);
template<> const GeomSymmetricTensor<3> GeomSymmetricTensor<3>::one = GeomSymmetricTensor<3>(1.0, 0.0, 0.0,
                                                                                             0.0, 1.0, 0.0,
                                                                                             0.0, 0.0, 1.0);
template<> const GeomSymmetricTensor<3>::size_type GeomSymmetricTensor<3>::mNumElements = 6;

template<> const double GeomSymmetricTensor<1>::onethird = 1.0/3.0;
template<> const double GeomSymmetricTensor<2>::onethird = 1.0/3.0;
template<> const double GeomSymmetricTensor<3>::onethird = 1.0/3.0;

template<> const double GeomSymmetricTensor<1>::sqrt3 = std::sqrt(3.0);
template<> const double GeomSymmetricTensor<2>::sqrt3 = std::sqrt(3.0);
template<> const double GeomSymmetricTensor<3>::sqrt3 = std::sqrt(3.0);

}

