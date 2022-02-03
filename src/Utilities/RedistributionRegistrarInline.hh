namespace Spheral {

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
inline
RedistributionRegistrar&
RedistributionRegistrar::
instance() {
  static RedistributionRegistrar theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// The non-const iterators.
//------------------------------------------------------------------------------
inline
RedistributionRegistrar::iterator
RedistributionRegistrar::
begin() {
   return mRedistributionNotificationHandles.begin();
}

inline
RedistributionRegistrar::iterator
RedistributionRegistrar::
end() {
   return mRedistributionNotificationHandles.end();
}

//------------------------------------------------------------------------------
// The const iterators.
//------------------------------------------------------------------------------
inline
RedistributionRegistrar::const_iterator
RedistributionRegistrar::
begin() const {
   return mRedistributionNotificationHandles.begin();
}

inline
RedistributionRegistrar::const_iterator
RedistributionRegistrar::
end() const {
   return mRedistributionNotificationHandles.end();
}

}
