//---------------------------------Spheral++----------------------------------//
// AxialSymmetryBoundary
// This boundary condition enforces (1D) axial symmetry on all fields by 
// creating (slave) ghost nodes in concentric cylindrical shells about a 1-D 
// line of internal nodes.  This boundary condition may only be used with 3-D 
// simulations.
//
// Created by JNJ.
//----------------------------------------------------------------------------//
#ifndef AxialSymmetryBoundary_HH
#define AxialSymmetryBoundary_HH

#ifndef __GCCXML__
#include <map>
#else
#include "fakestl.hh"
#endif

#include "Boundary.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace BoundarySpace {

class AxialSymmetryBoundary: public Boundary<Dim<3> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::ThirdRankTensor ThirdRankTensor;

  typedef Boundary<Dim<3> >::BoundaryNodes BoundaryNodes;

  // Constructors and destructors.
  AxialSymmetryBoundary(KernelSpace::TableKernel<Dim<3> >* kernel);
  virtual ~AxialSymmetryBoundary();

  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dim<3> >& nodeList);
  virtual void updateGhostNodes(NodeSpace::NodeList<Dim<3> >& nodeList);

  // Find the set of nodes in violation of this boundary in the given NodeList.
  // For planar boundaries this is any node that is "behind" the enter plane.
  virtual void setViolationNodes(NodeSpace::NodeList<Dim<3> >& nodeList);
  virtual void updateViolationNodes(NodeSpace::NodeList<Dim<3> >& nodeList);

  // Apply the boundary condition to ghost nodes.
  virtual void applyGhostBoundary(FieldSpace::Field<Dim<3> , int>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dim<3> , Scalar>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dim<3> , Vector>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dim<3> , Tensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dim<3> , SymTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dim<3>, ThirdRankTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dim<3>, std::vector<Scalar> >& field) const;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dim<3> , int>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dim<3> , Scalar>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dim<3> , Vector>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dim<3> , Tensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dim<3> , SymTensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dim<3>, ThirdRankTensor>& field) const;

private:

  KernelSpace::TableKernel<Dim<3> >* mKernel;
  std::map<NodeSpace::NodeList<Dim<3> >*, double> mR; // Ghost radii for nodelists.

  // Generic ghost-boundary applier.
  template <typename Value>
  void mApplyGhostBoundary(FieldSpace::Field<Dim<3>, Value>& field) const;
};

}
}

#else

namespace Spheral {
namespace BoundarySpace {
// Forward declaration.
class AxialSymmetryBoundary;
}
}

#endif
