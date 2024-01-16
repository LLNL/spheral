//---------------------------------Spheral++----------------------------------//
// A varient on the Monaghan-Gingold pairwise viscosity.
// This form sums the MG pairwise interactions for each node ahead of time,
// creating an overall viscous pressure for each node.  Therefore this is like 
// the vonNeuman viscosity using the MG method for determining velocity gradients.
// References: 
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// Created by JMO, Sun May 21 23:46:02 PDT 2000
//----------------------------------------------------------------------------//
#ifndef MonaghanGingoldSumViscosity_HH
#define MonaghanGingoldSumViscosity_HH

#include "MonaghanGingoldViscosity.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
  template<typename Dimension, typename DataType> class FieldList;
}

namespace Spheral {

template<typename Dimension>
class MonaghanGingoldSumViscosity: public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef NodeIDIterator<Dimension> IDIterator;
  typedef typename vector<Boundary<Dimension>*>::const_iterator ConstBoundaryIterator;

  // Constructors.
  MonaghanGingoldSumViscosity();
  MonaghanGingoldSumViscosity(Scalar Clinear, Scalar Cquadratic);

  // Destructor.
  ~MonaghanGingoldSumViscosity();

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual void initialize(const DataBase<Dimension>& dataBase,
                          ConstBoundaryIterator boundaryBegin,
                          ConstBoundaryIterator boundaryEnd,
                          const Scalar time,
                          const Scalar dt,
                          const TableKernel<Dimension>& W);

  // Method to calculate and return the viscous acceleration, work, and pressure,
  // all in one step (efficiency and all).
  virtual void viscousEffects(Vector& acceleration,
                              Scalar& work,
                              Scalar& pressure,
                              const IDIterator& nodeI,
                              const IDIterator& nodeJ,
                              const Vector& rij, 
                              const Vector& vi, const Vector& vj,
                              const Vector& etai, const Vector& etaj,
                              const Scalar ci, const Scalar cj,
                              const Scalar rhoi, const Scalar rhoj,
                              const Vector& gradW) const;

  // Access the viscous internal energy.
  const FieldList<Dimension, Scalar>& viscousInternalEnergyField() const;

  // Test if the ArtificialViscosity is valid.
  virtual bool valid() const;

private:
  //--------------------------- Private Interface ---------------------------//
  FieldList<Dimension, Scalar> mViscousEnergy;
};

}

#include "Field/NodeIDIterator.hh"

namespace Spheral {

using Spheral::NodeIDIterator;

//------------------------------------------------------------------------------
// Method to calculate and return the viscous acceleration, work, and pressure,
// all in one step (efficiency and all).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
MonaghanGingoldSumViscosity<Dimension>::
viscousEffects(Vector& acceleration,
               Scalar& work,
               Scalar& pressure,
               const IDIterator& nodeI,
               const IDIterator& nodeJ,
               const Vector& rij, 
               const Vector& vi, const Vector& vj,
               const Vector& etai, const Vector& etaj,
               const Scalar ci, const Scalar cj,
               const Scalar rhoi, const Scalar rhoj,
               const Vector& gradW) const {
  REQUIRE(valid());
  const Vector vij = vi - vj;
  const Scalar Qepsi = viscousInternalEnergyField()(nodeI);
  CHECK(Qepsi*rhoi >= 0.0);
  const Scalar QPi = Qepsi/rhoi;
  acceleration = QPi*gradW;
  work = QPi*vij.dot(gradW);
  pressure = rhoi*Qepsi;
}

//------------------------------------------------------------------------------
// Return the viscous energy field list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MonaghanGingoldSumViscosity<Dimension>::
viscousInternalEnergyField() const {
  return mViscousEnergy;
}

}

#endif
