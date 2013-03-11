//---------------------------------Spheral++----------------------------------//
// CylindricalBoundary -- Create a 3-D boundary around a set of nodes in the
// xy plane to emulate cylindrical (RZ) coordinates.
//
// Created by JMO, Tue Mar 29 13:48:03 PST 2005
//----------------------------------------------------------------------------//

#ifndef CylindricalBoundary_HH
#define CylindricalBoundary_HH

#include "Boundary.hh"
#include "Geometry/Dimension.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataOutput/registerWithRestart.hh"

namespace Spheral {
  namespace FileIOSpace {
    class FileIO;
  }
}

namespace Spheral {
namespace BoundarySpace {

class CylindricalBoundary: public Boundary<Dim<3> > {

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
  CylindricalBoundary(const DataBaseSpace::DataBase<Dim<3> >& dataBase);
  virtual ~CylindricalBoundary();

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

  // Find the effective reflection operator between the given xy
  // position and slaved ghost position.
  static Tensor reflectOperator(const Vector& r0, const Vector& r1);

  // Compute the target angular spacing for the given position and H.
  static double angularSpacing(const double ri,
                               const double hzi,
                               const double nodePerh,
                               const double kernelExtent);

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "CylindricalBoundary"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  FieldSpace::FieldList<Dim<3>, Scalar> mDeltaPhi;
  FieldSpace::FieldList<Dim<3>, Vector> mGhostPositions;

#ifndef __GCCXML__
  // The restart registration.
  DataOutput::RestartRegistrationType mRestart;
#endif
};

}
}

#ifndef __GCCXML__
#include "CylindricalBoundaryInline.hh"
#endif

#else

namespace Spheral {
  namespace BoundarySpace {
    // Forward declaration.
    class CylindricalBoundary;
  }
}

#endif
