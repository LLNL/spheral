//---------------------------------Spheral++----------------------------------//
// AxisBoundary -- A specialized Reflecting boundary condition for use along
// the axis in RZ calculations.
// Note: this boundary is automatically constructred by the RZ hydro objects,
//       so the user should not explicitly add this boundary.
//
// Created by JMO, Sun Aug 14 10:20:53 PDT 2016
//----------------------------------------------------------------------------//
#ifndef AxisBoundary_HH
#define AxisBoundary_HH

#include "Geometry/Dimension.hh"
#include "ReflectingBoundary.hh"

namespace Spheral {
namespace BoundarySpace {

class AxisBoundaryRZ: public ReflectingBoundary<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  AxisBoundaryRZ(const double etamin);
  virtual ~AxisBoundaryRZ();

  // Find the set of nodes in violation of this boundary in the given NodeList.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Access the fuzz from the axis we're using to enforce the BC.
  double etamin() const;
  void etamin(const double x);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "AxisBoundaryRZ"; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  double mEtaMin;
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    class AxisBoundaryRZ;
  }
}

#endif
