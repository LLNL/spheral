namespace Spheral {

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
inline
const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>&
SolidCRKSPHHydroBaseRZ::
deviatoricStressTT() const {
  return mDeviatoricStressTT;
}

inline
const FieldSpace::FieldList<Dim<2>, Dim<2>::Scalar>&
SolidCRKSPHHydroBaseRZ::
DdeviatoricStressTTDt() const {
  return mDdeviatoricStressTTDt;
}

}
