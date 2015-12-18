//---------------------------------Spheral++----------------------------------//
// FileIO -- Provide the generic interface to file objects.
//
// Created by J. Michael Owen, Wed Feb  7 22:59:14 PST 2001
//----------------------------------------------------------------------------//

#include <cstring>
#include "boost/algorithm/string.hpp"

#include "FileIO.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace FileIOSpace {

using namespace std;
using FieldSpace::Field;
using FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Empty Constructor.
//------------------------------------------------------------------------------
FileIO::FileIO():
  mFileName(""),
  mAccess(Undefined),
  mFileOpen(false)
{
}

//------------------------------------------------------------------------------
// Construct with the given file name and access type.
//------------------------------------------------------------------------------
FileIO::FileIO(const string filename, AccessType access):
  mFileName(filename),
  mAccess(access),
  mFileOpen(false)
{
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
FileIO::~FileIO()
{
}

//------------------------------------------------------------------------------
// Split the given path into its component strings.
//------------------------------------------------------------------------------
vector<string>
FileIO::splitPathComponents(const string path) const {

  vector<string> components;
  boost::split(components, path, boost::is_any_of("/"));

  // string currentComponent = "";
  // for (string::const_iterator itr = path.begin(); itr < path.end(); ++itr) {
  //   currentComponent += *itr;
  //   if (*itr == '/') {
  //     components.push_back(currentComponent);
  //     currentComponent = "";
  //   }
  // }
  // if (*(components.end() - 1) != currentComponent) components.push_back(currentComponent);

  return components;
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
// Provide access to the Scalar write method with float
//------------------------------------------------------------------------------
void
FileIO::write(const float& value, const string pathName) {
  write(double(value), pathName);
}

//------------------------------------------------------------------------------
// Provide access to the Scalar read method with float
//------------------------------------------------------------------------------
void
FileIO::read(float& value, const string pathName) const {
  double scalarValue;
  read(scalarValue, pathName);
  value = float(scalarValue);
}

//------------------------------------------------------------------------------
// Provide access to the vector<Scalar> write method with vector<float>
//------------------------------------------------------------------------------
void
FileIO::write(const vector<float>& value, const string pathName) {
  vector<double> scalarValue;
  scalarValue.reserve(value.size());
  for (vector<float>::const_iterator valItr = value.begin();
       valItr != value.end();
       ++valItr) scalarValue.push_back(*valItr);
  write(scalarValue, pathName);
}

//------------------------------------------------------------------------------
// Provide access to the vector<Scalar> read method with vector<float>
//------------------------------------------------------------------------------
void
FileIO::read(vector<float>& value, const string pathName) const {
  vector<double> scalarValue;
  read(scalarValue, pathName);
  value.resize(0);
  value.reserve(scalarValue.size());
  for (vector<double>::const_iterator valItr = scalarValue.begin();
       valItr != scalarValue.end();
       ++valItr) value.push_back(*valItr);
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

//   PyObject* result = PyObject_CallMethod(mPickleMod, "loads", "s", thpt);

//   cerr << "      result is ";
//   PyObject_Print(result, stderr, 0);
//   cerr << "   from   ";
//   cerr << "           "
//        << encodedThing 
//        << endl;

  // Deallocate stuff.
  Py_DECREF(pyEncodedThing);
  Py_DECREF(argsToPickle);

  // We're done.
  Py_INCREF(result);
  return result;
}
#endif

}

}
