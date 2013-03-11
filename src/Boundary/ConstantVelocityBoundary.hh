//---------------------------------Spheral++----------------------------------//
// ConstantVelocityBoundary -- A boundary condition to enforce a constant 
// velocity on a given set of nodes.
//
// This boundary is very specialized -- it explicitly works on only one 
// NodeList.
//
// Created by JMO, Mon Sep  9 23:39:39 PDT 2002
//----------------------------------------------------------------------------//
#ifndef ConstantVelocityBoundary_HH
#define ConstantVelocityBoundary_HH

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "Boundary.hh"
#include "DataOutput/registerWithRestart.hh"

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
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class ConstantVelocityBoundary: public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  ConstantVelocityBoundary(const NodeSpace::NodeList<Dimension>& nodeList,
                           const std::vector<int>& nodeIndicies);
  virtual ~ConstantVelocityBoundary();

  //**********************************************************************
  // All Boundary conditions must provide the following methods:
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Scalar> >& field) const;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;
  //**********************************************************************

  // Allow read only access to the node indicies and their forced velocities.
  const NodeSpace::NodeList<Dimension>& nodeList() const;
  std::vector<int> nodeIndicies() const;
  std::vector<Vector> velocityCondition() const;

  // Determine if the boundary is in a "valid", ready to use state.
  virtual bool valid() const;

  //******************************************************************************
  // Restart methods.
  virtual std::string label() const { return "ConstantVelocityBoundary"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //******************************************************************************

protected:
  //--------------------------- Protected Interface ---------------------------//

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  const NodeSpace::NodeList<Dimension>* mNodeListPtr;
  FieldSpace::Field<Dimension, int> mNodes;
  FieldSpace::Field<Dimension, Vector> mVelocity;

  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif
};

}
}

#ifndef __GCCXML__
#include "ConstantVelocityBoundaryInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    template<typename Dimension> class ConstantVelocityBoundary;
  }
}

#endif
