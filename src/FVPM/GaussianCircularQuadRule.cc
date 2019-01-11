//---------------------------------Spheral++----------------------------------//
// GaussianCircularQuadRule
//
// Created by JNJ, Sun Jul 11 12:50:51 PDT 2010
//----------------------------------------------------------------------------//

#include "GaussianCircularQuadRule.hh"
#include "gauss_legendre.hh"


namespace Spheral {

//-------------------------------------------------------------------
GaussianCircularQuadRule<Dim<1> >::
GaussianCircularQuadRule(const TableKernel<Dim<1> >& W, size_t N):
  CircularQuadRule<Dim<1> >(W),
  mX(N),
  mW(N),
  mA(-1.0),
  mB(1.0)
{
  // Generate the Gaussian quadrature table. Legendre polynomial 
  // zeros are generated to 1e-10 precision if they're not already in 
  // the table.
  gauss_legendre_tbl(static_cast<int>(N), &mX[0], &mW[0], 1e-10);
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
GaussianCircularQuadRule<Dim<1> >::
~GaussianCircularQuadRule()
{
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void 
GaussianCircularQuadRule<Dim<1> >::
setDomains(const Dim<1>::Vector& x1, 
           double r1,
           const Dim<1>::Vector& x2,
           double r2)
{
  // Find the interval [a, b].
  mB = x1.x() + r1;
  mA = mB - (x1.x() - r2);
  if (x1.x() > x2.x())
    swap(mA, mB);
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void 
GaussianCircularQuadRule<Dim<1> >::
getPointsAndWeights(vector<Dim<1>::Vector>& points,
                    vector<double>& weights) const
{
  REQUIRE(points.size() == weights.size());
  for (size_t q = 1; q < points.size(); ++q)
  {
    points[q].x(0.5 * (mB + mA) + 0.5 * (mB - mA) * mX[q]);
    weights[q] = 0.5 * (mB - mA) * mW[q];
  }
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
GaussianCircularQuadRule<Dim<2> >::
GaussianCircularQuadRule(const TableKernel<Dim<2> >& W):
  CircularQuadRule<Dim<2> >(W)
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
GaussianCircularQuadRule<Dim<2> >::
~GaussianCircularQuadRule()
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void 
GaussianCircularQuadRule<Dim<2> >::
setDomains(const Dim<2>::Vector& x1, 
           double r1,
           const Dim<2>::Vector& x2,
           double r2)
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void 
GaussianCircularQuadRule<Dim<2> >::
getPointsAndWeights(vector<Dim<2>::Vector>& points,
                    vector<double>& weights) const
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
GaussianCircularQuadRule<Dim<3> >::
GaussianCircularQuadRule(const TableKernel<Dim<3> >& W):
  CircularQuadRule<Dim<3> >(W)
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
GaussianCircularQuadRule<Dim<3> >::
~GaussianCircularQuadRule()
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void 
GaussianCircularQuadRule<Dim<3> >::
setDomains(const Dim<3>::Vector& x1, 
           double r1,
           const Dim<3>::Vector& x2,
           double r2)
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
void 
GaussianCircularQuadRule<Dim<3> >::
getPointsAndWeights(vector<Dim<3>::Vector>& points,
                    vector<double>& weights) const
{
  CHECK(false); // Nope!!
}
//-------------------------------------------------------------------

}
