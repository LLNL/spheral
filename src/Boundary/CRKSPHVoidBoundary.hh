//---------------------------------Spheral++----------------------------------//
// CRKSPHVoidBoundary -- Special boundary class to create void points off of
// the surface.
//
// Created by JMO, Fri Aug 18 10:31:50 PDT 2017
//
// Modified by:
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPHVoidBoundary__
#define __Spheral_CRKSPHVoidBoundary__

#include "Boundary.hh"

namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace MeshSpace {
    template<typename Dimension> class Mesh;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class CRKSPHVoidBoundary : public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  CRKSPHVoidBoundary(const FieldSpace::FieldList<Dimension, int>& surfacePoint,
                     const FieldSpace::FieldList<Dimension, std::vector<Vector>>& etaVoidPoints);
  virtual ~CRKSPHVoidBoundary();

  //**********************************************************************
  // All Boundary conditions must provide the following methods:
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList) override;

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the ghost node values in the given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Scalar> >& field) const override;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList) override;

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const override;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const override;
  //**********************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  const FieldSpace::FieldList<Dimension, int>& mSurfacePoint;
  const FieldSpace::FieldList<Dimension, std::vector<Vector>>& mEtaVoidPoints;
};

}
}

#else

namespace Spheral {
  namespace BoundarySpace {
    // Forward declaration.
    template<typename Dimension> class CRKSPHVoidBoundary;
  }
}

#endif
