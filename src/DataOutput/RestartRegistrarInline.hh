namespace Spheral {

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
inline
RestartRegistrar&
RestartRegistrar::
instance() {
  static RestartRegistrar theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// The non-const iterators.
//------------------------------------------------------------------------------
inline
RestartRegistrar::iterator
RestartRegistrar::
begin() {
   return mRestartHandles.begin();
}

inline
RestartRegistrar::iterator
RestartRegistrar::
end() {
   return mRestartHandles.end();
}

//------------------------------------------------------------------------------
// The const iterators.
//------------------------------------------------------------------------------
inline
RestartRegistrar::const_iterator
RestartRegistrar::
begin() const {
   return mRestartHandles.begin();
}

inline
RestartRegistrar::const_iterator
RestartRegistrar::
end() const {
   return mRestartHandles.end();
}

}
