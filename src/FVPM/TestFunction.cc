//---------------------------------Spheral++----------------------------------//

#include "TestFunction.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

//-------------------------------------------------------------------
template <typename Dimension>
TestFunction<Dimension>::
TestFunction(const TableKernel<Dimension>& W,
             const QuadRule<Dimension>& quadRule):
  mKernel(W),
  mQuadRule(quadRule)
{
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template <typename Dimension>
TestFunction<Dimension>::
~TestFunction()
{
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template <typename Dimension>
double
TestFunction<Dimension>::
operator()(const typename Dimension::Vector& xi,
           const typename Dimension::SymTensor& Hi,
           const vector<typename Dimension::Vector>& xjs,
           const vector<typename Dimension::SymTensor>& Hjs,
           const typename Dimension::Vector& x) const
{
  Vector etai = Hi.dot(x - xi);
  double etaiMag = etai.magnitude();
  double detHi = Hi.det();
  double Wi = mKernel.kernelValue(etai, detHi);

  double sumWj = Wi;
  for (size_t j = 0; j < xjs.size(); ++j)
  {
    Vector etaj = Hi.dot(xjs[j] - xi);
    double etajMag = etaj.magnitude();
    double detHj = Hjs[j].det();
    double Wj = mKernel.kernelValue(etajMag, detHj);
    sumWj += Wj;
  }

  return Wi/sumWj;
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template <typename Dimension>
typename Dimension::Vector
TestFunction<Dimension>::
interactionVector(const typename Dimension::Vector& xi,
                  const typename Dimension::SymTensor& Hi,
                  const typename Dimension::Vector& xj,
                  const typename Dimension::SymTensor& Hj,
                  const vector<typename Dimension::Vector>& xks,
                  const vector<typename Dimension::SymTensor>& Hks) const
{
  // Get the quadrature points and weights from the quadrature rule on the 
  // domain centered at xi and having the extent Hi.
  mQuadRule.setDomains(xi, Hi, xj, Hj);
  vector<Vector> points;
  vector<double> weights;
  mQuadRule.getPointsAndWeights(xi, Hi);

  // Now integrate.
  Vector betaij;
  for (size_t q = 0; q < points.size(); ++q)
  {
    Vector xq = points[q];

    // Compute Wi and its gradient at xq.
    Vector etai = Hi.dot(xi - xq);
    double etaiMag = etai.magnitude();
    double detHi = Hi.det();
    double Wi = mKernel.kernelValue(etai, detHi);
    Vector Hetai = Hi.dot(etai);
    Vector gradWi = mKernel.gradValue(etai, detHi) * Hetai;

    // Compute Wj and its gradient at xq.
    Vector etaj = Hi.dot(xj - xq);
    double etajMag = etaj.magnitude();
    double detHj = Hj.det();
    double Wj = mKernel.kernelValue(etaj, detHj);
    Vector Hetaj = Hj.dot(etaj);
    Vector gradWj = mKernel.gradValue(etaj, detHj) * Hetaj;

    // Now compute the sum of the kernels of the nodes {k}.
    double sumWk = 0.0;
    for (size_t k = 0; k < xks.size(); ++k)
    {
      Vector etak = Hks[k].dot(xks[k] - xq);
      double etakMag = etak.magnitude();
      double detHk = Hks[k].det();
      double Wk = mKernel.kernelValue(etakMag, detHk);
      sumWk += Wk;
    }

    // Contribute to the integral.
    double weight = weights[q];
    betaij += weight * (Wi * gradWj - Wj * gradWi) / (sumWk * sumWk);
  }

  return betaij;
}
//-------------------------------------------------------------------

}
