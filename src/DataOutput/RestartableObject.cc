//---------------------------------Spheral++----------------------------------//
// RestartableObject
// This is an object that handles registering it's descendents with the 
// RestartRegistrar.  The idea here is you just inherit from this bad boy, 
// and blamo you're restartable.
// This is intended solely for use making Python objects play in our restart
// setup -- please use the C++ centric methods in registerWithRestart.hh
// for C++ restarting.
//
// Created by JMO, Thu May 28 17:47:48 PDT 2009
//----------------------------------------------------------------------------//
#include "RestartableObject.hh"
#include "FileIO/FileIO.hh"

// These types are defined in our pybindgen generated code.
#ifndef _PyBindGenWrapperFlags_defined_
#define _PyBindGenWrapperFlags_defined_
typedef enum _PyBindGenWrapperFlags {
   PYBINDGEN_WRAPPER_FLAG_NONE = 0,
   PYBINDGEN_WRAPPER_FLAG_OBJECT_NOT_OWNED = (1<<0),
} PyBindGenWrapperFlags;
#endif

typedef struct {
    PyObject_HEAD
    Spheral::FileIO *obj;
    PyObject *inst_dict;
    PyBindGenWrapperFlags flags:8;
} PySpheralFileIO;

extern PyTypeObject PySpheralFileIO_Type;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
RestartableObject::
RestartableObject(PyObject* self,
                  const unsigned priority):
  mRestart(registerWithRestart(*this, priority)),
  mSelf(self) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
RestartableObject::
~RestartableObject() {
}

//------------------------------------------------------------------------------
// label
//------------------------------------------------------------------------------
std::string
RestartableObject::
label() const {
  PyObject* pylabel = PyObject_CallMethod(mSelf, (char*) "label", NULL);
  VERIFY2(pylabel != NULL,
          "RestartableObject::label ERROR: must provide label method.");
  VERIFY2(PyString_Check(pylabel),
          "RestartableObject::label ERROR: label method must return a string.");
  std::string result = PyString_AsString(pylabel);
  Py_DECREF(pylabel);
  return result;
}

//------------------------------------------------------------------------------
// dumpState
//------------------------------------------------------------------------------
void
RestartableObject::
dumpState(FileIO& file, const std::string pathName) const {
  PySpheralFileIO* py_file = PyObject_GC_New(PySpheralFileIO, &PySpheralFileIO_Type);
  py_file->inst_dict = NULL;
  py_file->obj = &file;
  PyObject* result = PyObject_CallMethod(mSelf, (char*) "dumpState", (char*) "Os", py_file, pathName.c_str());
  VERIFY2(result != NULL,
          "RestartableObject::dumpState ERROR encountered calling python dumpState method.");
}

//------------------------------------------------------------------------------
// restoreState
//------------------------------------------------------------------------------
void
RestartableObject::
restoreState(const FileIO& file, const std::string pathName) {
  PySpheralFileIO* py_file = PyObject_GC_New(PySpheralFileIO, &PySpheralFileIO_Type);
  py_file->inst_dict = NULL;
  py_file->obj = const_cast<FileIO*>(&file);
  PyObject* result = PyObject_CallMethod(mSelf, (char*) "restoreState", (char*) "Os", py_file, pathName.c_str());
  VERIFY2(result != NULL,
          "RestartableObject::restoreState ERROR encountered calling python restoreState method.");
}

}
