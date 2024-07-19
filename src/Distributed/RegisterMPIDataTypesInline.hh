namespace Spheral {

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
inline
RegisterMPIDataTypes&
RegisterMPIDataTypes::
instance() {
  static RegisterMPIDataTypes theInstance;
  return theInstance;
}

}
