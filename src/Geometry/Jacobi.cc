// An adaptation of the Numerical Recipes Jacobi routine to calculate the eigenvalues
// and eigenvectors of an NxN matrix.  Slightly rewritten here to use C++ constructs.

// ** This routine assumes that underflows are set to zero! **

#include <cmath>
#include <vector>
using namespace std;

#include "Utilities/DBC.hh"

namespace NumericalRecipes {

#define ROTATE(a, i, j, k, l) g = a[i][j]; h = a[k][l]; a[i][j] = g - s*(h + g*tau);\
  a[k][l] = h + s*(g - h*tau);

int jacobi(vector< vector<float> >& a,
           vector<float>& d,
           vector< vector<float> >& v) {

  int n = a.size();
  d.resize(n, 0.0);
  v.resize(n);
  for (int i = 0; i < n; ++i) {
    CHECK(a[i].size() == n);
    v[i].resize(n, 0.0);
  }

  int j, iq, ip, i;
  float tresh, theta, tau, t, sm, s, h, g, c;
  vector<float> b(n, 0.0);
  vector<float> z(n, 0.0);

  for (ip = 0; ip < n; ++ip) {
    for (iq = 0; iq < n; ++iq) v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < n; ++ip) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  int nrot = 0;
  for (i = 1; i <= 50; ++i) {
    sm = 0.0;
    for (ip = 0; ip < n-1; ++ip) {
      for (iq = ip + 1; iq < n; ++iq) {
        sm += fabs(a[ip][iq]);
      }
    }
    if (sm == 0.0) return nrot;

    if (i < 4) {
      tresh = 0.2*sm/(n*n);
    } else {
      tresh = 0.0;
    }

    for (ip = 0; ip < n-1; ++ip) {
      for (iq = ip + 1; iq < n; ++iq) {
        g = 100.0*fabs(a[ip][iq]);
        if (i > 4 && (float) (fabs(d[ip]) + g) == (float) fabs(d[ip]) &&
            (float) (fabs(d[iq]) + g) == (float) (fabs(d[iq]))) {
          a[ip][iq] = 0.0;
        } else if (fabs(a[ip][iq]) > tresh) {
          h = d[iq] - d[ip];
          if ((float) (fabs(h) + g) == (float) fabs(h)) {
            t = (a[ip][iq])/h;
          } else {
            theta = 0.5*h/(a[ip][iq]);
            t = 1.0/(fabs(theta) + sqrt(1.0 + theta*theta));
            if (theta < 0.0) t = -t;
          }
          c = 1.0/sqrt(1 + t*t);
          s = t*c;
          tau = s/(1.0 + c);
          h = t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq] = 0.0;
          for (j = 0; j < ip-1; ++j) {
            ROTATE(a, j, ip, j, iq);
          }
          for (j = ip + 1; j <= iq - 1; ++j) {
            ROTATE(a, ip, j, j, iq);
          }
          for (j = iq + 1; j < n; ++j) {
            ROTATE(a, ip, j, iq, j);
          }
          for (j = 0; j < n; ++j) {
            ROTATE(v, j, ip, j, iq);
          }
          ++nrot;
        }
      }
    }

    for (ip = 0; ip < n; ++ip) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  throw "Too many iterations in routine jacobi.";
}
}
