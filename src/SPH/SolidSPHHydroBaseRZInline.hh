namespace Spheral {

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
inline
const FieldList<Dim<2>, Dim<2>::Scalar>&
SolidSPHHydroBaseRZ::
deviatoricStressTT() const {
  return mDeviatoricStressTT;
}

inline
const FieldList<Dim<2>, Dim<2>::Scalar>&
SolidSPHHydroBaseRZ::
DdeviatoricStressTTDt() const {
  return mDdeviatoricStressTTDt;
}

}
