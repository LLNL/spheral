//---------------------------------Spheral++----------------------------------//
// PeriodicBoundary -- Apply a Periodic boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Apr 19 14:00:50 PDT 2000
//----------------------------------------------------------------------------//

#ifndef PeriodicBoundary_HH
#define PeriodicBoundary_HH

#include "PlanarBoundary.hh"

namespace Spheral {
namespace BoundarySpace {

class PeriodicPlanarBoundary;

template<typename Dimension>
class PeriodicBoundary: public PlanarBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  PeriodicBoundary();
  PeriodicBoundary(const GeomPlane<Dimension>& plane1,
                   const GeomPlane<Dimension>& plane2);
  virtual ~PeriodicBoundary();

  // Override the default methods for setting ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Override the default methods for setting violation nodes.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Override the methods for setting the enter and exit planes.
  virtual const GeomPlane<Dimension>& enterPlane() const;
  virtual void setEnterPlane(const GeomPlane<Dimension>& enterPlane);

  virtual const GeomPlane<Dimension>& exitPlane() const;
  virtual void setExitPlane(const GeomPlane<Dimension>& exitPlane);

  // Override the culling of ghost nodes.
  virtual void cullGhostNodes(const FieldSpace::FieldList<Dimension, int>& flagSet,
                              FieldSpace::FieldList<Dimension, int>& old2newIndexMap,
                              std::vector<int>& numNodesRemoved);

  // Apply the boundary condition to the ghost nodes in the given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Scalar> >& field) const;

  // Enforce the boundary condition on the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;

  // Override the base reset method.
  virtual void reset(const DataBaseSpace::DataBase<Dimension>& dataBase);

  // Report the number of ghost nodes in this boundary.
  virtual int numGhostNodes() const;

  virtual std::string label() const { return "PeriodicBoundary"; }

private:
  //--------------------------- Private Interface ---------------------------//
  // The PeriodicBounary actually is made up of two nested class boundary 
  // conditions, which we hide from the user for simplicity of interface.
//   class PeriodicPlanarBoundary;
#include "PeriodicPlanarBoundary.hh"

  PeriodicPlanarBoundary mPlane1Boundary;
  PeriodicPlanarBoundary mPlane2Boundary;
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    template<typename Dimension> class PeriodicBoundary;
  }
}

#endif
