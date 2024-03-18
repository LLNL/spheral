//---------------------------------Spheral++----------------------------------//
// An experimental hour glass control algorithm for SPH, based on the ASPH
// second moment ideas.
//
// Created by JMO, Sun Jan 15 21:19:53 PST 2006
//----------------------------------------------------------------------------//
#include "SecondMomentHourglassControl.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Inline functions specialized to Dimension to calculate the anti-hourglassing
// acceleration.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
pairwiseAntihourglassing(const typename Dimension::Scalar& /*rji*/,
                         const typename Dimension::Scalar& /*Wi*/,
                         const typename Dimension::Scalar& /*Wj*/,
                         const typename Dimension::Scalar& /*gradWi*/,
                         const typename Dimension::Scalar& /*gradWj*/,
                         const double /*dt2*/) {
  VERIFY(false);
  return 0.0;
}

template<>
inline
Dim<1>::Scalar
pairwiseAntihourglassing<Dim<1> >(const Dim<1>::Scalar& rji,
                                  const Dim<1>::Scalar& Wi,
                                  const Dim<1>::Scalar& Wj,
                                  const Dim<1>::Scalar& gradWi,
                                  const Dim<1>::Scalar& gradWj,
                                  const double dt2) {
  CONTRACT_VAR(rji);
  REQUIRE(rji >= 0.0);
  REQUIRE(Wi + Wj >= 0.0);
  REQUIRE(dt2 >= 0.0);
  const double tiny = 1.0e-30;
  const double denom = Wi*gradWi - Wj*gradWj;
  return (Wi*Wi - Wj*Wj)/(denom*dt2 + sgn(denom)*tiny);
}

template<>
inline
Dim<2>::Scalar
pairwiseAntihourglassing<Dim<2> >(const Dim<2>::Scalar& rji,
                                  const Dim<2>::Scalar& Wi,
                                  const Dim<2>::Scalar& Wj,
                                  const Dim<2>::Scalar& gradWi,
                                  const Dim<2>::Scalar& gradWj,
                                  const double dt2) {
  REQUIRE(rji >= 0.0);
  REQUIRE(Wi + Wj >= 0.0);
  REQUIRE(dt2 >= 0.0);
  const double tiny = 1.0e-30;
  const double thpt = Wi*Wi - Wj*Wj;
  const double denom = Wi*gradWi - Wj*gradWj - thpt/(rji + tiny);
  return thpt/(denom*dt2 + sgn(denom)*tiny);
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SecondMomentHourglassControl<Dimension>::
SecondMomentHourglassControl(const TableKernel<Dimension>& W,
                             const double multiplier,
                             const double maxAccelerationFactor):
  Physics<Dimension>(),
  mW(W),
  mMultiplier(multiplier),
  mMaxAccelerationFactor(maxAccelerationFactor),
  mAcceleration(FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SecondMomentHourglassControl<Dimension>::
~SecondMomentHourglassControl() {
}

//------------------------------------------------------------------------------
// Determine the principle derivatives for the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SecondMomentHourglassControl<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  //const double tiny = 1.0e-30;
  const double dt2 = dt*dt;

  // Get the state fields.
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, SymTensor> Hfield = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, Scalar());

  // Prepare to fill in the diagnostic acceleration field.
  mAcceleration = dataBase.newFluidFieldList(Vector::zero, "anti-hourglass acceleration");

  // Get the connectivity map.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();

  // Iterate over the NodeLists.
  for (auto iNodeList = 0u; iNodeList != nodeLists.size(); ++iNodeList) {
    const FluidNodeList<Dimension>* nodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(nodeLists[iNodeList]);
    CHECK(nodeListPtr != 0);
    const Field<Dimension, Vector>& r = **position.fieldForNodeList(*nodeListPtr);
    //const Field<Dimension, Vector>& v = **velocity.fieldForNodeList(*nodeListPtr);
    const Field<Dimension, SymTensor>& H = **Hfield.fieldForNodeList(*nodeListPtr);
    Field<Dimension, Vector>& accel = **DvDt.fieldForNodeList(*nodeListPtr);
    //Field<Dimension, Scalar>& work = **DepsDt.fieldForNodeList(*nodeListPtr);

    // The diagnostic acceleration.
    Field<Dimension, Vector>& diagnostic = **mAcceleration.fieldForNodeList(*nodeListPtr);

    // Iterate over the nodes in this NodeList.
    for (auto i = 0u; i != nodeListPtr->numInternalNodes(); ++i) {

      // State for node i.
      const Vector& ri = r(i);
      //const Vector& vi = v(i);
      const SymTensor& Hi = H(i);
      const Scalar Hdeti = Hi.Determinant();

      // Find the neighbors for this node.
      // Note that we only want neighbors from the same NodeList...
      const vector<int>& neighbors = connectivityMap.connectivityForNode(nodeListPtr, i)[iNodeList];

      // Iterate over the neighbors, and build up the vote for the hourglass motion
      Vector hg;
      for (vector<int>::const_iterator jItr = neighbors.begin();
           jItr != neighbors.end();
           ++jItr) {
        const int j = *jItr;
        if ((int)i != j) {

          // State for node j.
          const Vector& rj = r(j);
          //const Vector& vj = v(j);
          const SymTensor& Hj = H(j);
          const Scalar Hdetj = Hj.Determinant();

          // Compute the acceleration from this pair.
          const Vector rij = ri - rj;
          const Vector rijUnit = rij.unitVector();
          const Scalar rijMag = rij.magnitude();

          // Compute the rotation necessary to align with the radial vector
          // between these nodes.
          const Tensor Ri = rotationMatrix(rijUnit);

          // Kernel estimate from i.
          const Scalar etai = (Hi*rij).magnitude();
          const Scalar Wi = mW(etai, Hdeti);
          const Scalar gradWi = mW.grad(etai, Hdeti);

          // Kernel estimate from j.
          const Scalar etaj = (Hj*rij).magnitude();
          const Scalar Wj = mW(etaj, Hdetj);
          const Scalar gradWj = mW.grad(etaj, Hdetj);

          // Acceleration.
          const Scalar ai = 0.5*pairwiseAntihourglassing<Dimension>(rijMag,
                                                                    Wi,
                                                                    Wj,
                                                                    gradWj,
                                                                    gradWi,
                                                                    dt2);

          // We only use this acceleration if it will heat the system -- dissipitive 
          // processes only please!
          const Vector hgij = ai*rijUnit;
//           if (vi.dot(hgij) < 0.0) 
            hg += hgij;

        }
      }

      // Apply the correction.
      hg *= mMultiplier;
      const Vector DvDti = min(mMaxAccelerationFactor*(accel(i).magnitude()), hg.magnitude()) * hg.unitVector();
//       const double fac = min(1.0, mMaxAccelerationFactor*vi.magnitude() / (0.5*hg.magnitude()*dt + tiny));
//       const Vector DvDti = fac*hg;
      //const Scalar worki = -(vi.dot(DvDti));
//       if (worki > 0.0) {
        diagnostic(i) += DvDti;
        accel(i) += DvDti;
//         work(i) += worki;
//       }

    }
  }

}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename SecondMomentHourglassControl<Dimension>::TimeStepType
SecondMomentHourglassControl<Dimension>::
dt(const DataBase<Dimension>&, 
   const State<Dimension>&,
   const StateDerivatives<Dimension>&,
   const typename Dimension::Scalar) const {
  return TimeStepType(FLT_MAX, "No vote.");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SecondMomentHourglassControl<Dimension>::
registerState(DataBase<Dimension>&,
              State<Dimension>&) {
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SecondMomentHourglassControl<Dimension>::
registerDerivatives(DataBase<Dimension>&,
                    StateDerivatives<Dimension>&) {
}

}

