#include <string>

// BPL includes
#include "boost/python.hpp"
#include "boost/python/str.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;

//------------------------------------------------------------------------------
// Print a string.
//------------------------------------------------------------------------------
const char*
printSTDString(const std::string& self) {
  return self.c_str();
}

//------------------------------------------------------------------------------
// Some useful methods for strings.
//------------------------------------------------------------------------------
boost::python::str
operator+(const std::string& lhs,
	  const boost::python::str& rhs) {
  return boost::python::str(lhs).join(rhs);
}

boost::python::str
operator+(const boost::python::str& lhs,
	  const std::string& rhs) {
  return lhs.join(boost::python::str(rhs));
}

//------------------------------------------------------------------------------
// Automatically convert C++ string <--> Python str.
// This is directly cribbed from Boost.Python documentation.
//------------------------------------------------------------------------------
// template <>
// object AutoConverter<std::string>::toObject( std::string const& s )
// {
//   return str( s.cstr() );
// }

// template <>
// void* AutoConverter<std::string>::convertible( PyObject* obj_ptr )
// {
//   if ( !PyString_Check(obj_ptr) ) return 0;
//   return obj_ptr;
// }

// template <>
// void AutoConverter<std::string>::fromPython(PyObject* obj_ptr,
// 					    void* memblock)
// {
//   const char* value = PyString_AsString( obj_ptr );
//   if (value == 0) {
//     throw_error_already_set();
//   }
//   new ( memblock ) std::string( value );
// }


// struct string_to_python_str {
//   static PyObject* convert(string const& s) {
//     return boost::python::incref(boost::python::object(s).ptr());
//   }
// };

// struct string_from_python_str {
//   string_from_python_str() {
//     boost::python::converter::registry::push_back(&convertible,
//                                                   &construct,
//                                                   boost::python::type_id<string>());
//   }

//   static void* convertible(PyObject* obj_ptr) {
//     if (!PyString_Check(obj_ptr)) return 0;
//     return obj_ptr;
//   }

//   static void construct(PyObject* obj_ptr,
//                         boost::python::converter::rvalue_from_python_stage1_data* data) {
//     const char* value = PyString_AsString(obj_ptr);
//     if (value == 0) boost::python::throw_error_already_set();
//     void* storage = ((boost::python::converter::rvalue_from_python_storage<string>*)
//                      data)->storage.bytes;
//     new (storage) string(value);
//     data->convertible = storage;
//   }
// };

// std::size_t size(string const& s) { return s.size(); }

//------------------------------------------------------------------------------
// Wrap the std::string.
//------------------------------------------------------------------------------
void
wrapSTDString() {

//   boost::python::to_python_converter<string, string_to_python_str>();
//   string_from_python_str();

  class_<std::string>("string", init<>())

    // Constructors.
    .def(init<char*>())

    .def("__len__", &std::string::size)

    // String representations.
    .def("__str__", printSTDString)
    .def("__repr__", printSTDString)

    .def(self + self)
//     .def( self + other<boost::python::str> )

    ;

//   AutoConverter<std::string>();
}


}
