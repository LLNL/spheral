namespace Spheral{

// Accessor Fns
template<typename Dimension>
inline
void
CullenDehnenViscosity<Dimension>::
alphMax(Scalar val)
{
    malphMax = val;
}

template<typename Dimension>
inline
void
CullenDehnenViscosity<Dimension>::
alphMin(Scalar val)
{
    malphMin = val;
}

template<typename Dimension>
inline
void
CullenDehnenViscosity<Dimension>::
betaE(Scalar val)
{
    mbetaE = val;
}

template<typename Dimension>
inline
void
CullenDehnenViscosity<Dimension>::
betaD(Scalar val)
{
    mbetaD = val;
}

template<typename Dimension>
inline
void
CullenDehnenViscosity<Dimension>::
betaC(Scalar val)
{
    mbetaC = val;
}

template<typename Dimension>
inline
void
CullenDehnenViscosity<Dimension>::
fKern(Scalar val)
{
    mfKern = val;
}

template<typename Dimension>
inline
void
CullenDehnenViscosity<Dimension>::
boolHopkins(bool val)
{
    mboolHopkins = val;
}

//------------------------------------------------------------------------------
// Access the main kernel
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const TableKernel<Dimension>&
CullenDehnenViscosity<Dimension>::kernel() const {
  return mKernel;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
CullenDehnenViscosity<Dimension>::PrevDvDt() const {
   return mPrevDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::PrevDivV() const {
   return mPrevDivV;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::PrevDivV2() const {
   return mPrevDivV2;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::CullAlpha() const {
   return mCullAlpha;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::CullAlpha2() const {
   return mCullAlpha2;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::R() const {
   return mR;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::vsig() const {
   return mVsig;
}

template<typename Dimension>
inline
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
alphMax() const{return malphMax;}

template<typename Dimension>
inline
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
alphMin() const{return malphMin;}

template<typename Dimension>
inline
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
betaE() const{return mbetaE;}

template<typename Dimension>
inline
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
betaD() const{return mbetaD;}

template<typename Dimension>
inline
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
betaC() const{return mbetaC;}

template<typename Dimension>
inline
typename Dimension::Scalar
CullenDehnenViscosity<Dimension>::
fKern() const{return mfKern;}

template<typename Dimension>
inline
bool 
CullenDehnenViscosity<Dimension>::
boolHopkins() const{return mboolHopkins;}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::
DalphaDt() const{ return mDalphaDt;}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
CullenDehnenViscosity<Dimension>::
alphaLocal() const{ return mAlphaLocal;}

}