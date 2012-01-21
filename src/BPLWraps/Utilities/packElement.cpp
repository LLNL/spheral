#include <boost/python.hpp>
#include <string>
#include <vector>
#include <algorithm>

#include "Field/Field.hh"
#include "Utilities/packElement.hh"

using namespace boost::python;
using namespace std;

namespace Spheral {

//------------------------------------------------------------------------------
// Wrapper around our packElement method more convenient for exposing to Python.
//------------------------------------------------------------------------------
template<typename Element>
inline
string
convertElementToString(const Element& x) {
  vector<char> buffer;
  packElement(x, buffer);
  return string(&(*buffer.begin()), buffer.size());
}

//------------------------------------------------------------------------------
// Wrapper around our unpackElement method more convenient for exposing to
//  Python.
//------------------------------------------------------------------------------
template<typename Element>
inline
Element
convertStringToElement(const string& x) {
  vector<char> buffer;
  copy(x.begin(), x.end(), back_inserter(buffer));
  Element result;
  vector<char>::const_iterator itr = buffer.begin();
  unpackElement(result, itr, buffer.end());
  return result;
}

//------------------------------------------------------------------------------
// Wrap our useful little methods to converting various element types to string
// byte representations.
//------------------------------------------------------------------------------
void wrapPackElement() {
  def("packElementInt", &convertElementToString<int>);
  def("packElementUL", &convertElementToString<unsigned long>);
  def("packElementULL", &convertElementToString<unsigned long long>);
  def("packElementDouble", &convertElementToString<double>);

  def("unpackElementInt", &convertStringToElement<int>);
  def("unpackElementUL", &convertStringToElement<unsigned long>);
  def("unpackElementULL", &convertStringToElement<unsigned long long>);
  def("unpackElementDouble", &convertStringToElement<double>);
}

}
