namespace Spheral {
namespace SolidSPHSpace {

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
inline
const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>&
SolidSPHHydroBaseRZ::
deviatoricStressTT() const {
  return mDeviatoricStressTT;
}

inline
const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>&
SolidSPHHydroBaseRZ::
DdeviatoricStressTTDt() const {
  return mDdeviatoricStressTTDt;
}

}
}
