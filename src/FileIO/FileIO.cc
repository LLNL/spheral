//---------------------------------Spheral++----------------------------------//
// FileIO -- Provide the generic interface to file objects.
//
// Created by J. Michael Owen, Wed Feb  7 22:59:14 PST 2001
//----------------------------------------------------------------------------//
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

#ifndef CXXONLY
namespace py = pybind11;
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Empty Constructor.
//------------------------------------------------------------------------------
FileIO::FileIO():
  mFileName(""),
  mAccess(AccessType::Undefined),
  mFileOpen(false) {
}

//------------------------------------------------------------------------------
// Construct with the given file name and access type.
//------------------------------------------------------------------------------
FileIO::FileIO(const string filename, AccessType access):
  mFileName(filename),
  mAccess(access),
  mFileOpen(false) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
FileIO::~FileIO() {
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
// Write polytopes
//------------------------------------------------------------------------------
void
FileIO::write(const Dim<1>::FacetedVolume& value, const string pathName) {
  std::vector<char> buf;
  packElement(value, buf);
  std::string bufstr(buf.begin(), buf.end());
  this->write(bufstr, pathName);
}

void
FileIO::write(const Dim<2>::FacetedVolume& value, const string pathName) {
  std::vector<char> buf;
  packElement(value, buf);
  std::string bufstr(buf.begin(), buf.end());
  this->write(bufstr, pathName);
}

void
FileIO::write(const Dim<3>::FacetedVolume& value, const string pathName) {
  std::vector<char> buf;
  packElement(value, buf);
  std::string bufstr(buf.begin(), buf.end());
  this->write(bufstr, pathName);
}

//------------------------------------------------------------------------------
// Read polytopes
//------------------------------------------------------------------------------
void
FileIO::read(Dim<1>::FacetedVolume& value, const string pathName) const {
  std::string bufstr;
  this->read(bufstr, pathName);
  const std::vector<char> buf(bufstr.begin(), bufstr.end());
  auto itr = buf.begin();
  unpackElement(value, itr, buf.end());
  ENSURE(itr == buf.end());
}

void
FileIO::read(Dim<2>::FacetedVolume& value, const string pathName) const {
  std::string bufstr;
  this->read(bufstr, pathName);
  const std::vector<char> buf(bufstr.begin(), bufstr.end());
  auto itr = buf.begin();
  unpackElement(value, itr, buf.end());
  ENSURE(itr == buf.end());
}

void
FileIO::read(Dim<3>::FacetedVolume& value, const string pathName) const {
  std::string bufstr;
  this->read(bufstr, pathName);
  const std::vector<char> buf(bufstr.begin(), bufstr.end());
  auto itr = buf.begin();
  unpackElement(value, itr, buf.end());
  ENSURE(itr == buf.end());
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
// Python objects (handle with pickle)
//------------------------------------------------------------------------------
void
FileIO::write_object(py::object thing, const std::string& pathName) {
  auto pickle = py::module_::import("pickle");
  auto dumps = pickle.attr("dumps");
  auto stuff = dumps(thing);
  //std::cerr << "-- writeObject string size, bytes size: " << std::string(stuff).size() << " " << py::len(stuff) << std::endl;
  this->write_bytes(stuff, pathName);
}

py::object
FileIO::read_object(const std::string& pathName) const {
  auto pickle = py::module_::import("pickle");
  auto loads = pickle.attr("loads");
  auto stuff = this->read_bytes(pathName);
  return loads(stuff);
}

//------------------------------------------------------------------------------
// Python bytes
//------------------------------------------------------------------------------
void
FileIO::write_bytes(py::bytes stuff, const std::string& pathName) {
  this->write(std::string(stuff), pathName);
}

py::bytes
FileIO::read_bytes(const std::string& pathName) const {
  std::string stuff;
  this->read(stuff, pathName);
  return py::bytes(stuff);
}
#endif

}
