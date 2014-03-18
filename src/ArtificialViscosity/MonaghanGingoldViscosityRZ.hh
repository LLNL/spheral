//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
// This version specialized for cylindrical (RZ) coordinates.
// References: 
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// Created by JMO, Mon Nov 20 15:50:29 PST 2006
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
  MonaghanGingoldViscosityRZ(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  virtual ~MonaghanGingoldViscosityRZ();

  // Method to apply the viscous acceleration, work, and pressure, to the derivatives
  // all in one step (efficiency and all).
  virtual void viscousEffects(const DataBaseSpace::DataBase<Dimension>& dataBase,
                              const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                              const State<Dimension>& state,
                              StateDerivatives<Dimension>& derivatives) const;

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
