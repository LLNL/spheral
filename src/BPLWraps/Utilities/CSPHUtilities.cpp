#include <boost/python.hpp>
#include "Utilities/CSPHUtilities.hh"
#include "Geometry/Dimension.hh"

using namespace boost::python;

namespace Spheral {

//------------------------------------------------------------------------------
// Make a version of CSPHKernelAndGradient that returns the result as a tuple.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
tuple
CSPHKernelAndGradientPY(const KernelSpace::TableKernel<Dimension>& W,
                        const typename Dimension::Vector& rij,
                        const typename Dimension::Vector& etaj,
                        const typename Dimension::SymTensor& Hj,
                        const typename Dimension::Scalar& Hdetj,
                        const typename Dimension::Scalar& Ai,
                        const typename Dimension::Vector& Bi,
                        const typename Dimension::Vector& gradAi,
                        const typename Dimension::Tensor& gradBi) {
  typename Dimension::Scalar WCSPH;
  typename Dimension::Scalar gradWSPH;
  typename Dimension::Vector gradWCSPH;
  CSPHKernelAndGradient(W, rij, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, WCSPH, gradWSPH, gradWCSPH);
  return make_tuple(WCSPH, gradWSPH, gradWCSPH);
}

//------------------------------------------------------------------------------
// Expose the methods for computing the CSPH corrections.
//------------------------------------------------------------------------------
void wrapCSPHUtilities() {
  def("CSPHKernel", &CSPHKernel<Dim<1> >, "Compute the CSPH corrected kernel estimate.");
  def("CSPHKernel", &CSPHKernel<Dim<2> >, "Compute the CSPH corrected kernel estimate.");
  def("CSPHKernel", &CSPHKernel<Dim<3> >, "Compute the CSPH corrected kernel estimate.");

  def("CSPHKernelAndGradient", &CSPHKernelAndGradientPY<Dim<1> >, "Compute the CSPH corrected kernel and gradient estimates.");
  def("CSPHKernelAndGradient", &CSPHKernelAndGradientPY<Dim<2> >, "Compute the CSPH corrected kernel and gradient estimates.");
  def("CSPHKernelAndGradient", &CSPHKernelAndGradientPY<Dim<3> >, "Compute the CSPH corrected kernel and gradient estimates.");
}

}
