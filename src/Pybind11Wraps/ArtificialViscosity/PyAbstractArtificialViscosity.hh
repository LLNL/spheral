//------------------------------------------------------------------------------
// Trampoline classes for abstract ArtificialViscosity interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyAbstractArtificialViscosity__
#define __Spheral_PyAbstractArtificialViscosity__

#include <utility>

#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "DataBase/DataBase.hh"

#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;

namespace Spheral {
namespace ArtificialViscositySpace {

//------------------------------------------------------------------------------
// PyAbstractArtificialViscosity
//------------------------------------------------------------------------------
template<typename Dimension, class ArtificialViscosityBase>
class PyAbstractArtificialViscosity: public ArtificialViscosityBase {
public:
  using ArtificialViscosityBase::ArtificialViscosityBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename ArtificialViscosityBase::ConstBoundaryIterator ConstBoundaryIterator;

  virtual void initialize(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
                          const Scalar time,
                          const Scalar dt,
                          const KernelSpace::TableKernel<Dimension>& W) override {
    PYBIND11_OVERLOAD(void,                    // Return type
                      ArtificialViscosityBase, // Parent class
                      initialize,              // name of method
                      dataBase, state, derivs, boundaryBegin, boundaryEnd, time, dt, W      // arguments
      );
  }

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
    PYBIND11_OVERLOAD_PURE(ReturnType,                                        // Return type
                           ArtificialViscosityBase,                           // Parent class
                           Piij,                                              // name of method
                           nodeListi, i, nodeListj, j, xi, etai, vi, rhoi, csi, Hi, xj, etaj, vj, rhoj, csj, Hj       // arguments
      );
  }

};

}
}

#endif
