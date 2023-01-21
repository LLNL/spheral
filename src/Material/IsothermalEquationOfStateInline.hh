namespace Spheral {

//------------------------------------------------------------------------------
// Access the polytropic constant.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
K() const {
  return mK;
}

//------------------------------------------------------------------------------
// Access the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
molecularWeight() const {
  return mMolecularWeight;
}

}
