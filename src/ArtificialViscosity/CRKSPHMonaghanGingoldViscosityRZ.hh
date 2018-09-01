//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is specialized for use with CRKSPH.
//
// This form specialized for use with the area-weighted RZ formalism.
//
// Created by JMO, Sun May 22 10:45:30 PDT 2016
//----------------------------------------------------------------------------//
#ifndef CRKSPHMonaghanGingoldViscosityRZ_HH
#define CRKSPHMonaghanGingoldViscosityRZ_HH

#include "CRKSPHMonaghanGingoldViscosity.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

class CRKSPHMonaghanGingoldViscosityRZ: public CRKSPHMonaghanGingoldViscosity<Dim<2> > {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Dimension::FourthRankTensor FourthRankTensor;
  typedef Dimension::FifthRankTensor FifthRankTensor;
  typedef ArtificialViscosity<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CRKSPHMonaghanGingoldViscosityRZ(const Scalar Clinear,
                                   const Scalar Cquadratic,
                                   const bool linearInExpansion,
                                   const bool quadraticInExpansion,
                                   const Scalar etaCritFrac,
                                   const Scalar etaFoldFrac);

  // Destructor.
  virtual ~CRKSPHMonaghanGingoldViscosityRZ();

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
  virtual std::string label() const { return "CRKSPHMonaghanGingoldViscosityRZ"; }

private:
  //--------------------------- Private Interface ---------------------------//
  CRKSPHMonaghanGingoldViscosityRZ();
  CRKSPHMonaghanGingoldViscosityRZ(const CRKSPHMonaghanGingoldViscosityRZ&);
  CRKSPHMonaghanGingoldViscosityRZ& operator=(const CRKSPHMonaghanGingoldViscosityRZ&) const;
};

}

#else

namespace Spheral {
  // Forward declaration.
  class CRKSPHMonaghanGingoldViscosityRZ;
}

#endif
