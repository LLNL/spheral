//---------------------------------Spheral++----------------------------------//
// SphericalBoundary -- Create a 3-D boundary around a line of nodes down
// the x-axis appropriate for enforcing a spherically symmetric system.
//
// Created by JMO, Tue Mar 15 21:39:43 PST 2005
//----------------------------------------------------------------------------//

#ifndef SphericalBoundary_HH
#define SphericalBoundary_HH

#include "Boundary.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Geometry/Dimension.hh"
#include "DataOutput/registerWithRestart.hh"

namespace FileIOSpace {
  class FileIO;
}

namespace Spheral {
namespace BoundarySpace {

class SphericalBoundary: public Boundary<Dim<3> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<3> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;
  typedef Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Boundary<Dimension>::BoundaryNodes BoundaryNodes;

  // Constructors and destructors.
  SphericalBoundary(const DataBaseSpace::DataBase<Dim<3> >& dataBase);
  virtual ~SphericalBoundary();

  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the ghost node values in the given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Scalar> >& field) const;

  // Find the set of nodes in violation of this boundary in the given NodeList.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;

  // Find the effective reflection operator between the given x-axis
  // position and slaved ghost position.
  Tensor reflectOperator(const Vector& r0, const Vector& r1) const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "SphericalBoundary"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  FieldSpace::FieldList<Dim<3>, std::vector<Vector> > mGhostPositions;

#ifndef __GCCXML__
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif
};

}
}

#ifndef __GCCXML__
#include "SphericalBoundaryInline.hh"
#endif

#else

namespace Spheral {
  namespace BoundarySpace {
    // Forward declaration.
    class SphericalBoundary;
  }
}

#endif
