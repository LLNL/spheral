//---------------------------------Spheral++----------------------------------//
// SphericalOriginBoundary -- A specialized Reflecting boundary condition for
// use at the origin in Spherical coordinate calculations.
// Note: this boundary is automatically constructred by the Spherical hydro objects,
//       so the user should not explicitly add this boundary.
//
// Created by JMO, Fri Feb 25 15:03:39 PST 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalOriginBoundary__
#define __Spheral_SphericalOriginBoundary__

#include "Geometry/Dimension.hh"
#include "ReflectingBoundary.hh"

namespace Spheral {

class SphericalOriginBoundary: public ReflectingBoundary<Dim<1>> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<1> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  SphericalOriginBoundary();
  virtual ~SphericalOriginBoundary();

  // Override the ghost node methods to not make ghost nodes
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override;

  // Find the set of nodes in violation of this boundary in the given NodeList.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SphericalOriginBoundary"; }
  //****************************************************************************

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
};

}

#else

// Forward declaration.
namespace Spheral {
  class SphericalOriginBoundary;
}

#endif
