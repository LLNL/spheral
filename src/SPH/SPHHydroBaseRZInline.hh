namespace Spheral {
namespace SPHSpace {

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
inline
const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>&
SPHHydroBaseRZ::
massDensityRZ() const {
  return mMassDensityRZ;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>&
SPHHydroBaseRZ::
f1() const {
  return mf1;
}

template<typename Dimension>
inline
const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>&
SPHHydroBaseRZ::
f2() const {
  return mf2;
}

}
}
