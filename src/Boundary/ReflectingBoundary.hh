//---------------------------------Spheral++----------------------------------//
// ReflectingBoundary -- Apply a Reflecting boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//

#ifndef ReflectingBoundary_HH
#define ReflectingBoundary_HH

#include "Boundary.hh"
#include "PlanarBoundary.hh"

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class ReflectingBoundary: public PlanarBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  ReflectingBoundary();
  ReflectingBoundary(const GeomPlane<Dimension>& plane);
  virtual ~ReflectingBoundary();

  // Apply the boundary condition to the ghost values of given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;

  // Apply the boundary condition to face centered fields on a tessellation.
  virtual void enforceBoundary(std::vector<int>& faceField, const MeshSpace::Mesh<Dimension>& mesh) const;
  virtual void enforceBoundary(std::vector<Scalar>& faceField, const MeshSpace::Mesh<Dimension>& mesh) const;
  virtual void enforceBoundary(std::vector<Vector>& faceField, const MeshSpace::Mesh<Dimension>& mesh) const;
  virtual void enforceBoundary(std::vector<Tensor>& faceField, const MeshSpace::Mesh<Dimension>& mesh) const;
  virtual void enforceBoundary(std::vector<SymTensor>& faceField, const MeshSpace::Mesh<Dimension>& mesh) const;
  virtual void enforceBoundary(std::vector<ThirdRankTensor>& faceField, const MeshSpace::Mesh<Dimension>& mesh) const;

  // Fill in faces on this boundary with effective opposite face values.
  virtual void swapFaceValues(FieldSpace::Field<Dimension, std::vector<Scalar> >& field,
                              const MeshSpace::Mesh<Dimension>& mesh) const;
  virtual void swapFaceValues(FieldSpace::Field<Dimension, std::vector<Vector> >& field,
                              const MeshSpace::Mesh<Dimension>& mesh) const;

  // Allow read only access to the reflection operator.
  const Tensor& reflectOperator() const;

  // Valid test.
  virtual bool valid() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "ReflectingBoundary"; }
  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  Tensor mReflectOperator;
};

}
}

#ifndef __GCCXML__
#include "ReflectingBoundaryInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    template<typename Dimension> class ReflectingBoundary;
  }
}

#endif
