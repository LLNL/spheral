//------------------------------------------------------------------------------
// Trampoline classes for virtual ArtificialViscosity interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyArtificialViscosity__
#define __Spheral_PyArtificialViscosity__

#include "PyAbstractArtificialViscosity.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;

namespace Spheral {
namespace ArtificialViscositySpace {

//------------------------------------------------------------------------------
// PyArtificialViscosity
//------------------------------------------------------------------------------
template<typename Dimension, class ArtificialViscosityBase>
class PyArtificialViscosity: public ArtificialViscosityBase {
public:
  using ArtificialViscosityBase::ArtificialViscosityBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  virtual std::pair<Tensor, Tensor> Piij(const unsigned nodeListi, const unsigned i, 
                                         const unsigned nodeListj, const unsigned j,
                                         const Vector& xi,
                                         const Vector& etai,
                                         const Vector& vi,
                                         const Scalar rhoi,
                                         const Scalar csi,
                                         const SymTensor& Hi,
                                         const Vector& xj,
                                         const Vector& etaj,
                                         const Vector& vj,
                                         const Scalar rhoj,
                                         const Scalar csj,
                                         const SymTensor& Hj) const override {
    typedef std::pair<Tensor, Tensor> ReturnType;
    PYBIND11_OVERLOAD(ReturnType,                                        // Return type
                      ArtificialViscosityBase,                           // Parent class
                      Piij,                                              // name of method
                      nodeListi, i, nodeListj, j, xi, etai, vi, rhoi, csi, Hi, xj, etaj, vj, rhoj, csj, Hj       // arguments
      );
  }

};

}
}

#endif
