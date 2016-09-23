// An adaptation of the Numerical Recipes Jacobi routine to calculate the eigenvalues
// and eigenvectors of an NxN matrix.  Slightly rewritten here to use C++ constructs.

namespace NumericalRecipes {
  int jacobi(vector< vector<float> >& a,
             vector<float>& d,
             vector< vector<float> >& v);
}
