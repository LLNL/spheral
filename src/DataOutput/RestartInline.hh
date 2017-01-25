namespace Spheral {
namespace DataOutput {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Object>
Restart<Object>::
Restart(Object& object):
  mObjectPtr(&object) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Object>
Restart<Object>::
~Restart() {
}

//------------------------------------------------------------------------------
// Dump the objects state.
//------------------------------------------------------------------------------
template<typename Object>
inline
void
Restart<Object>::
dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const {
  mObjectPtr->dumpState(file, pathName);
}

//------------------------------------------------------------------------------
// Restore the objects state.
//------------------------------------------------------------------------------
template<typename Object>
inline
void
Restart<Object>::
restoreState(const FileIOSpace::FileIO& file, const std::string& pathName) const {
  mObjectPtr->restoreState(file, pathName);
}

//------------------------------------------------------------------------------
// Get the label for the object.
//------------------------------------------------------------------------------
template<typename Object>
inline
std::string
Restart<Object>::
label() const {
  return mObjectPtr->label();
}

}
}
