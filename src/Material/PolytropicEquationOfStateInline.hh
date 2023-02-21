namespace Spheral {

//------------------------------------------------------------------------------
// Access the polytropic constant.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
polytropicConstant() const {
  return mPolytropicConstant;
}

//------------------------------------------------------------------------------
// Access the polytropic index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
polytropicIndex() const {
  return mPolytropicIndex;
}

//------------------------------------------------------------------------------
// Access gamma ( (n+1)/n ).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
gamma() const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Access the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
molecularWeight() const {
  return mMolecularWeight;
}

}
