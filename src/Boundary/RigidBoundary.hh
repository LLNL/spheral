//---------------------------------Spheral++----------------------------------//
// RigidBoundary -- This is identical to the ReflectingBoundary except that 
// all vectors other than x and v are copied, not reflected, across the 
// boundary.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//

#ifndef RigidBoundary_HH
#define RigidBoundary_HH

#include "Boundary.hh"
#include "PlanarBoundary.hh"

namespace Spheral {

template<typename Dimension>
class RigidBoundary: public PlanarBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  RigidBoundary();
  RigidBoundary(const GeomPlane<Dimension>& plane);
  virtual ~RigidBoundary();

  // Apply the boundary condition to the ghost values of given Field.
  virtual void applyGhostBoundary(Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, int>& field) const;
  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const;

  // Allow read only access to the reflection operator.
  const Tensor& reflectOperator() const;

  // Valid test.
  virtual bool valid() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const { return "RigidBoundary"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  Tensor mReflectOperator;
};

}

#include "RigidBoundaryInline.hh"

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class RigidBoundary;
}

#endif
