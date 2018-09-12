//---------------------------------Spheral++----------------------------------//
// GradPressureViscosity.  A viscosity based on pairwise estimates of the
// pressure gradient.
//
// Created by JMO, Fri Dec 13 15:37:42 PST 2002
//----------------------------------------------------------------------------//
#include "GradPressureViscosity.hh"
#include "Kernel/TableKernel.hh"

#include "Utilities/DBC.hh"

#include <algorithm>

namespace Spheral {


//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GradPressureViscosity<Dimension>::
GradPressureViscosity():
  ArtificialViscosity<Dimension>(),
  mKernelPtr(0),
  mWeight(FieldStorageType::ReferenceFields) {
}

//------------------------------------------------------------------------------
// Construct with a single scalar, representing the pressure coefficient.
//------------------------------------------------------------------------------
template<typename Dimension>
GradPressureViscosity<Dimension>::
GradPressureViscosity(Scalar Clinear, Scalar Cquadratic):
  ArtificialViscosity<Dimension>(Clinear, Cquadratic),
  mKernelPtr(0),
  mWeight(FieldStorageType::ReferenceFields) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GradPressureViscosity<Dimension>::
~GradPressureViscosity() {
}

//------------------------------------------------------------------------------
// Intialization method, called at the beginning of a cycle.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GradPressureViscosity<Dimension>::
initialize(const DataBase<Dimension>& dataBase,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryBegin,
           typename ArtificialViscosity<Dimension>::ConstBoundaryIterator boundaryEnd,
           const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const TableKernel<Dimension>& W) {

  // Invoke the base class initialization.
  ArtificialViscosity<Dimension>::initialize(dataBase, 
					     boundaryBegin,
					     boundaryEnd,
					     time,
					     dt,
					     W);

  // Maintain a pointer to the kernel being used.
  mKernelPtr = &W;
  mWeight = dataBase.fluidWeight();
}

//------------------------------------------------------------------------------
// Calculate the viscous acceleration, work, and equivalent pressure.
//------------------------------------------------------------------------------
// template<typename Dimension>
// void
// GradPressureViscosity<Dimension>::
// viscousEffects(typename Dimension::Vector& acceleration,
//                typename Dimension::Scalar& work,
//                typename Dimension::Scalar& pressure,
//                const NodeIDIterator<Dimension>& nodeI,
//                const NodeIDIterator<Dimension>& nodeJ,
//                const typename Dimension::Vector& rij, 
//                const typename Dimension::Vector& vi,
//                const typename Dimension::Vector& vj,
//                const typename Dimension::Vector& etai,
//                const typename Dimension::Vector& etaj,
//                const typename Dimension::Scalar ci,
//                const typename Dimension::Scalar cj,
//                const typename Dimension::Scalar Pi,
//                const typename Dimension::Scalar Pj,
//                const typename Dimension::Scalar rhoi,
//                const typename Dimension::Scalar rhoj,
//                const typename Dimension::Vector& gradW) const {


//   REQUIRE(negligibleSoundSpeed() > 0.0);
//   REQUIRE(epsilon2() > 0.0);
//   REQUIRE(rhoi > 0.0);
//   REQUIRE(mKernelPtr != 0);

//   const Vector rijUnit = rij.unitVector();
//   const Scalar hinverse = (mHfield(nodeI)*rijUnit).magnitude();
//   CHECK(hinverse > 0.0);

//   acceleration = Cl()*hinverse/rhoi*(Pi - Pj)/(etai.magnitude() + epsilon2())*
//     0.25*(mWeight(nodeI) + mWeight(nodeJ))*((*mKernelPtr)(etai, mHfield(nodeI)) +
//                                             (*mKernelPtr)(etaj, mHfield(nodeJ)))*
//     etai.unitVector();
//   work = -(vi.dot(acceleration));
// //   if (work > 0.0) {
//     pressure = Cl()*Pi;
// //   } else {
// //     acceleration.Zero();
// //     work = 0.0;
// //     pressure = 0.0;
// //   }
// }

template<typename Dimension>
void
GradPressureViscosity<Dimension>::
viscousEffects(typename Dimension::Vector& acceleration,
               typename Dimension::Scalar& work,
               typename Dimension::Scalar& pressure,
               const NodeIteratorBase<Dimension>& nodeI,
               const NodeIteratorBase<Dimension>& nodeJ,
               const typename Dimension::Vector& rij, 
               const typename Dimension::Vector& rijUnit, 
               const typename Dimension::Vector& vi,
               const typename Dimension::Vector& vj,
               const typename Dimension::Vector& etai,
               const typename Dimension::Vector& etaj,
               const typename Dimension::Scalar ci,
               const typename Dimension::Scalar cj,
               const typename Dimension::Scalar Pi,
               const typename Dimension::Scalar Pj,
               const typename Dimension::Scalar rhoi,
               const typename Dimension::Scalar rhoj,
               const typename Dimension::Scalar hi,
               const typename Dimension::Scalar hj,
               const typename Dimension::Vector& gradW) const {


  REQUIRE(this->negligibleSoundSpeed() > 0.0);
  REQUIRE(this->epsilon2() > 0.0);

  const Vector vij = vi - vj;

//   const Scalar gradP = std::abs(Pi - Pj)/(etai.magnitude() + epsilon2());
//     std::min(1.0, vij.magnitude()/(ci + negligibleSoundSpeed()));
//   const Scalar csi = std::max(ci, negligibleSoundSpeed());
//   const Scalar csj = std::max(cj, negligibleSoundSpeed());
//   const Scalar safei = Pi + epsilon2()*rhoi*csi*csi;
//   const Scalar safej = Pj + epsilon2()*rhoj*csj*csj;
//   const Scalar safe = safei + safej;

  const double Cl = this->Cl();

  const Scalar gradPi = std::abs(Pi - Pj)/(etai.magnitude() + this->epsilon2());
  const Scalar Qi = Cl*gradPi;
  const Scalar QPii = Qi/(rhoi*rhoi);

  const Scalar gradPj = std::abs(Pi - Pj)/(etaj.magnitude() + this->epsilon2());
  const Scalar Qj = Cl*gradPj;
  const Scalar QPij = Qj/(rhoj*rhoj);

  const Scalar Qij = 0.5*(Qi + Qj);
  const Scalar QPiij = 0.5*(QPii + QPij);

//   const Scalar QPi = Qi/(rhoi*rhoi) + Qi/(rhoj*rhoj);
//   const Scalar QPi = Qi/pow(min(rhoi, rhoj), 2);

  work = QPiij*vij.dot(gradW);
  if (work > 0.0) {
    acceleration = QPiij*gradW;
    pressure = Qij;
  } else {
    acceleration.Zero();
    work = 0.0;
    pressure = 0.0;
  }

}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class GradPressureViscosity< Dim<1> >;
  template class GradPressureViscosity< Dim<2> >;
  template class GradPressureViscosity< Dim<3> >;
}
