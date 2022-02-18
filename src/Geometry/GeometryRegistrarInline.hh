namespace Spheral {
  
//------------------------------------------------------------------------------
// Instance
//------------------------------------------------------------------------------
inline
GeometryRegistrar&
GeometryRegistrar::instance() {
  static GeometryRegistrar theInstance;
  return theInstance;
}

}
