//---------------------------------Spheral++----------------------------------//
// FileIO -- Provide the generic interface to file objects.
//
// Created by J. Michael Owen, Wed Feb  7 22:59:14 PST 2001
//----------------------------------------------------------------------------//
#ifndef CXXONLY
#include "Python.h"
#endif

#include "FileIO.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include "boost/algorithm/string.hpp"

#include <cstring>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Empty Constructor.
//------------------------------------------------------------------------------
FileIO::FileIO():
  mFileName(""),
  mAccess(AccessType::Undefined),
  mFileOpen(false)
#ifndef CXXONLY
  ,
  mPickleMod(NULL),
  mPickleDumps(NULL),
  mPickleLoads(NULL) {

  // Import pickle.
  PyObject* modName = Py_BuildValue("s", "cPickle");
  mPickleMod = PyImport_Import(modName);
  mPickleDumps = PyObject_GetAttrString(mPickleMod, "dumps");
  mPickleLoads = PyObject_GetAttrString(mPickleMod, "loads");
  Py_DECREF(modName);
#else
  {
#endif
}

//------------------------------------------------------------------------------
// Construct with the given file name and access type.
//------------------------------------------------------------------------------
FileIO::FileIO(const string filename, AccessType access):
  mFileName(filename),
  mAccess(access),
  mFileOpen(false)
#ifndef CXXONLY
  ,
  mPickleMod(NULL),
  mPickleDumps(NULL),
  mPickleLoads(NULL) {
  
  // Import pickle.
  PyObject* modName = Py_BuildValue("s", "cPickle");
  mPickleMod = PyImport_Import(modName);
  mPickleDumps = PyObject_GetAttrString(mPickleMod, "dumps");
  mPickleLoads = PyObject_GetAttrString(mPickleMod, "loads");
  Py_DECREF(modName);
#else
  {
#endif
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
FileIO::~FileIO() {
#ifndef CXXONLY
  Py_DECREF(mPickleMod);
  Py_DECREF(mPickleDumps);
  Py_DECREF(mPickleLoads);
#endif
}

//------------------------------------------------------------------------------
// Split the given path into its component strings.
//------------------------------------------------------------------------------
vector<string>
FileIO::splitPathComponents(const string path) const {
  vector<string> components;
  boost::split(components, path, boost::is_any_of("/"));
  return components;
}

// Reverse operation
std::string
FileIO::joinPathComponents(const std::vector<std::string>& components) const {
  string result = "";
  for (const auto s: components) result += "/" + s;
  return result;
}

//------------------------------------------------------------------------------
// Write GeomPlanes.
//------------------------------------------------------------------------------
void
FileIO::write(const GeomPlane<Dim<1> >& value, const string pathName) {
  write(value.point(), pathName + "/point");
  write(value.normal(), pathName + "/normal");
}

void
FileIO::write(const GeomPlane<Dim<2> >& value, const string pathName) {
  write(value.point(), pathName + "/point");
  write(value.normal(), pathName + "/normal");
}

void
FileIO::write(const GeomPlane<Dim<3> >& value, const string pathName) {
  write(value.point(), pathName + "/point");
  write(value.normal(), pathName + "/normal");
}

//------------------------------------------------------------------------------
// Read GeomPlanes.
//------------------------------------------------------------------------------
void
FileIO::read(GeomPlane<Dim<1> >& value, const string pathName) const {
  Dim<1>::Vector point, normal;
  read(point, pathName + "/point");
  read(normal, pathName + "/normal");
  value.point(point);
  value.normal(normal);
}

void
FileIO::read(GeomPlane<Dim<2> >& value, const string pathName) const {
  Dim<2>::Vector point, normal;
  read(point, pathName + "/point");
  read(normal, pathName + "/normal");
  value.point(point);
  value.normal(normal);
}

void
FileIO::read(GeomPlane<Dim<3> >& value, const string pathName) const {
  Dim<3>::Vector point, normal;
  read(point, pathName + "/point");
  read(normal, pathName + "/normal");
  value.point(point);
  value.normal(normal);
}

//------------------------------------------------------------------------------
// Write uniform_random
//------------------------------------------------------------------------------
void
FileIO::write(const uniform_random& value, const string pathName) {
  write((unsigned) value.seed(), pathName + "/seed");
  write((unsigned) value.numCalls(), pathName + "/numCalls");
  write(value.min(), pathName + "/min");
  write(value.max(), pathName + "/max");
}

//------------------------------------------------------------------------------
// Read uniform_random
//------------------------------------------------------------------------------
void
FileIO::read(uniform_random& value, const string pathName) const {
  unsigned seed, numCalls;
  double a, b;
  read(seed, pathName + "/seed");
  read(numCalls, pathName + "/numCalls");
  read(a, pathName + "/min");
  read(b, pathName + "/max");
  value.seed(seed);
  value.range(a, b);
  value.advance(numCalls);
}

//------------------------------------------------------------------------------
// Provide access to the string write method with char*
//------------------------------------------------------------------------------
void
FileIO::write(const char* value, const string pathName) {
  write(string(value, strlen(value)), pathName);
}

//------------------------------------------------------------------------------
// Provide access to the string read method with char*
//------------------------------------------------------------------------------
void
FileIO::read(char* value, const string pathName) const {
  string strValue;
  read(strValue, pathName);
  strcpy(value, strValue.c_str());
}

//------------------------------------------------------------------------------
// Return the group (directory) component of the given path.
//------------------------------------------------------------------------------
string
FileIO::groupName(const string pathName) const {

  const vector<string> components = splitPathComponents(pathName);
  CHECK(components.size() > 0);

  string groupName = "";
  for (vector<string>::const_iterator nameItr = components.begin();
       nameItr < components.end() - 1;
       ++nameItr) groupName += *nameItr;

  return groupName;
}

//------------------------------------------------------------------------------
// Return the variable name component of the given path.
//------------------------------------------------------------------------------
string
FileIO::variableName(const string pathName) const {
  const vector<string> components = splitPathComponents(pathName);
  CHECK(components.size() > 0);
  return components[components.size() - 1];
}

//------------------------------------------------------------------------------
// Return the current file name.
//------------------------------------------------------------------------------
const string&
FileIO::fileName() const {
  return mFileName;
}

//------------------------------------------------------------------------------
// Return the access being specified.
//------------------------------------------------------------------------------
AccessType
FileIO::access() const {
  return mAccess;
}

//------------------------------------------------------------------------------
// Return the flag indicating whether the file is open or not.
//------------------------------------------------------------------------------
bool
FileIO::fileOpen() const {
  return mFileOpen;
}

#ifndef CXXONLY
//------------------------------------------------------------------------------
// Write out a python object by pickling it.
//------------------------------------------------------------------------------
void
FileIO::
writeObject(PyObject* thing, PyObject* pathObj) {

  // Extract the path.
  const string path(PyString_AsString(pathObj));

  // Pickle our object.
  PyObject* thingArgs = Py_BuildValue("(O)", thing);
  PyObject* pickledThing = PyObject_CallObject(mPickleDumps, thingArgs);
  if (pickledThing == NULL) {
    PyErr_SetString(PyExc_ValueError, "Unable to pickle object");
  }

  // Extract the string representation of the object.
  CHECK(PyString_Check(pickledThing));
  string result(PyString_AsString(pickledThing));

  // Replace the \n to <<n>> to survive writing to the file.
  boost::replace_all(result, "\n", "<<n>>");
  CHECK(result.find("\n") == string::npos);

  // Now we can finally write the sucker out.
  this->write(result, string(path));

  // Deallocate stuff.
  Py_DECREF(thingArgs);
  Py_DECREF(pickledThing);
}

//------------------------------------------------------------------------------
// Read in and decode a pickled python object.
//------------------------------------------------------------------------------
PyObject*
FileIO::
readObject(PyObject* pathObj) const {

  // Extract the path.
  const string path(PyString_AsString(pathObj));

  // Read in the pickled string representation.
  string encodedThing;
  this->read(encodedThing, path);

  // Convert the \n's back to the real thing.
  boost::replace_all(encodedThing, "<<n>>", "\n");
  CHECK(encodedThing.find("<<n>>") == string::npos);

  // Turn the string into a python object.
  const char* thpt = encodedThing.c_str();
  PyObject* pyEncodedThing = Py_BuildValue("s", thpt);

  // Unpickle our object.
  PyObject* argsToPickle = Py_BuildValue("(O)", pyEncodedThing);
  PyObject* result = PyEval_CallObject(mPickleLoads, argsToPickle);
  if (result == NULL) {
    PyErr_SetString(PyExc_ValueError, "Unable to unpickle string");
    return NULL;
  }

  // cerr << "      result is ";
  // PyObject_Print(result, stderr, 0);
  // cerr << "   from   ";
  // cerr << "           "
  //      << encodedThing 
  //      << endl;

  // Deallocate stuff.
  Py_DECREF(pyEncodedThing);
  Py_DECREF(argsToPickle);

  // We're done.
  Py_INCREF(result);
  return result;
}
#endif

}
