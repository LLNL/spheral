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
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Constructors and destructors.
  PeriodicPlanarBoundary();
  PeriodicPlanarBoundary(const GeomPlane<Dimension>& plane1,
                         const GeomPlane<Dimension>& plane2);
  virtual ~PeriodicPlanarBoundary();

  // Apply the boundary condition to the ghost nodes in the given Field.
  virtual void applyGhostBoundary(Field<Dimension, int>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Enforce the boundary condition on the violation nodes in the given Field.
  virtual void enforceBoundary(Field<Dimension, int>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Valid test.
  virtual bool valid() const override;

private:
  //--------------------------- Private Interface ---------------------------//
};

#endif
