//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
// References: 
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// This form specialized for use with the area-weighted RZ formalism.
// 
// Created by JMO, Sat May 21 16:15:44 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_MonaghanGingoldViscosityRZ_hh__
#define __Spheral_MonaghanGingoldViscosityRZ_hh__

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {
namespace ArtificialViscositySpace {

class MonaghanGingoldViscosityRZ: public MonaghanGingoldViscosity<Dim<2> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  // Constructors.
  MonaghanGingoldViscosityRZ(const Scalar Clinear,
                             const Scalar Cquadratic,
                             const bool linearInExpansion,
                             const bool quadraticInExpansion);

  // Destructor.
  virtual ~MonaghanGingoldViscosityRZ();

  // The required method to compute the artificial viscous P/rho^2.
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
                                         const SymTensor& Hj) const;

  // Restart methods.
  virtual std::string label() const { return "MonaghanGingoldViscosityRZ"; }

private:
  //--------------------------- Private Interface ---------------------------//
  MonaghanGingoldViscosityRZ();
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace ArtificialViscositySpace {
    class MonaghanGingoldViscosityRZ;
  }
}

#endif
