//---------------------------------Spheral++----------------------------------//
// PeriodicPlanarBoundary -- The nested class member of Periodic Boundary that
// does the actual work.
//
// Created by JMO, Wed Apr 19 15:00:50 PDT 2000
//----------------------------------------------------------------------------//

#ifndef PeriodicPlanarBoundary_HH
#define PeriodicPlanarBoundary_HH

class PeriodicPlanarBoundary: public PlanarBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  PeriodicPlanarBoundary();
  PeriodicPlanarBoundary(const GeomPlane<Dimension>& plane1,
                         const GeomPlane<Dimension>& plane2);
  virtual ~PeriodicPlanarBoundary();

  // Apply the boundary condition to the ghost nodes in the given Field.
  virtual void applyGhostBoundary(Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const;

  // Enforce the boundary condition on the violation nodes in the given Field.
  virtual void enforceBoundary(Field<Dimension, int>& field) const;
  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const;

  // Valid test.
  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
};

#endif
