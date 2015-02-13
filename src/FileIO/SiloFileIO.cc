//---------------------------------Spheral++----------------------------------//
// SiloFileIO -- Provide the interface to silo file objects.
//
// Created by JMO, Sat Feb  7 23:06:03 PST 2015
//----------------------------------------------------------------------------//
#include <algorithm>
#include "boost/algorithm/string/replace.hpp"

#include "SiloFileIO.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace FileIOSpace {

using namespace std;
using FieldSpace::Field;

namespace {
//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
string SILO_mangle(const string& x) {
  string result = x;
  boost::replace_all(result, " ", "_");
  boost::replace_all(result, ".", "_p_");
  boost::replace_all(result, ",", "_c_");
  return result;
}
}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
SiloFileIO::SiloFileIO():
  FileIO(),
  mFilePtr(0) {
}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
SiloFileIO::
SiloFileIO(const string fileName, AccessType access):
  FileIO(fileName, access),
  mFilePtr(0) {
  open(fileName, access);
  ENSURE(mFileOpen && mFilePtr != 0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
SiloFileIO::~SiloFileIO() {
  close();
}

//------------------------------------------------------------------------------
// Open a SiloFile file with the specified access.
//------------------------------------------------------------------------------
void
SiloFileIO::open(const string fileName, AccessType access) {
  VERIFY2(mFilePtr == 0 and mFileOpen == false,
          "ERROR: attempt to reopen SiloFileIO object.");

  string fullFileName = fileName;
  if (fullFileName.find(".silo") == string::npos) {
    fullFileName += ".silo";
  }

  if (access == Read) {
    mFilePtr = DBOpen(fullFileName.c_str(), DB_HDF5, DB_READ);
  } else {
    mFilePtr = DBCreate(fullFileName.c_str(), DB_CLOBBER, DB_LOCAL, "Spheral++ restart file.", DB_HDF5);
  }
  VERIFY2(mFilePtr != 0, "SiloFileIO ERROR: unable to open " << fullFileName);
  mFileOpen = true;
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
void
SiloFileIO::close() {
  if (mFilePtr != 0) {
    VERIFY2(DBClose(mFilePtr) == 0,
            "SiloFileIO ERROR: unable to close file.");
    mFilePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const unsigned value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const int value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const bool value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  int ivalue = value ? 1 : 0;
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &ivalue, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const double value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const string value, const string pathName) {
  const char* cvalue = value.c_str();
  int dims[1] = {strlen(cvalue)};
  this->write(dims[0], pathName + "/size");
  if (dims[0] > 0) {
    const string varname = this->setDir(pathName + "/value");
    VERIFY2(DBWrite(mFilePtr, varname.c_str(), cvalue, dims, 1, DB_CHAR) == 0,
            "SiloFileIO ERROR: unable to write variable " << pathName);
  }
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Vector& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Tensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::SymTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::ThirdRankTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Vector& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {2};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Tensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {4};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::SymTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {3};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::ThirdRankTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {8};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Vector& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {3};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Tensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {9};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::SymTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {6};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::ThirdRankTensor& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {27};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Read an unsigned from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(unsigned& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(int& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), (void*) &value) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(bool& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  int ivalue;
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &ivalue) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
  value = (ivalue == 1 ? true : false);
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(double& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(string& value, const string pathName) const {
  int valsize;
  this->read(valsize, pathName + "/size");
  if (valsize == 0) {
    value = "";
  } else {
    const string varname = this->setDir(pathName + "/value");
    char* cvalue = (char*) DBGetVar(mFilePtr, varname.c_str());
    VERIFY2(cvalue != NULL,
            "SiloFileIO ERROR: unable to read variable " << pathName);
    value = string(cvalue);
    CHECK2(value.size() == valsize, value << " " << value.size() << " " << valsize);
    free(cvalue);
    // char cvalue[valsize + 1];  // Do we need to allow space for trailing null?
    // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), cvalue) == 0,
    //         "SiloFileIO ERROR: failed to read variable " << pathName);
    // value = string(cvalue);
  }
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Vector& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Tensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::SymTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::ThirdRankTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Vector& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Tensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::SymTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::ThirdRankTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Vector& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Tensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::SymTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::ThirdRankTensor& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<int>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {value.size()};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<double> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<double>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  int dims[1] = {value.size()};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<string> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<string>& value, const string pathName) {
  const unsigned n = value.size();
  vector<int> dim_stuff(n);
  string stuff;
  for (unsigned i = 0; i != n; ++i) {
    dim_stuff[i] = value[i].size();
    stuff += value[i];
  }
  this->write(dim_stuff, pathName + "/dim_stuff");
  this->write(stuff, pathName + "/stuff");

  // const string varname = this->setDir(pathName);
  // int dims[n];
  // char* cvalue[n];
  // for (unsigned i = 0; i != n; ++i) {
  //   dims[i] = value[i].size() + 1;
  //   CHECK2(strlen(value[i].c_str()) == dims[i]-1, i << " " << value[i] << " " << strlen(value[i].c_str()) << " " << dims[i]);
  //   cvalue[i] = new char[dims[i]];
  //   strcpy(cvalue[i], value[i].c_str());
  // }
  // VERIFY2(DBWrite(mFilePtr, varname.c_str(), cvalue, dims, 1, DB_CHAR) == 0,
  //         "SiloFileIO ERROR: unable to write variable " << pathName);
  // for (unsigned i = 0; i != n; ++i) {
  //   delete cvalue[i];
  // }
}

//------------------------------------------------------------------------------
// Write a vector<Dim<1>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<1>::Vector>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {n};
  vector<double> cvalue(n);
  for (unsigned i = 0; i != value.size(); ++i) cvalue[i] = value[i].x();
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<1>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<1>::Tensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {n};
  vector<double> cvalue(n);
  for (unsigned i = 0; i != value.size(); ++i) cvalue[i] = value[i].xx();
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<1>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<1>::SymTensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {n};
  vector<double> cvalue(n);
  for (unsigned i = 0; i != value.size(); ++i) cvalue[i] = value[i].xx();
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<1>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<1>::ThirdRankTensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {n};
  vector<double> cvalue(n);
  for (unsigned i = 0; i != value.size(); ++i) cvalue[i] = value[i](0,0,0);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<2>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<2>::Vector>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {2*n};
  vector<double> cvalue(2*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[2*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<2>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<2>::Tensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {4*n};
  vector<double> cvalue(4*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[4*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<2>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<2>::SymTensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {3*n};
  vector<double> cvalue(3*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[3*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<2>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<2>::ThirdRankTensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {8*n};
  vector<double> cvalue(8*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[8*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<3>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<3>::Vector>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {3*n};
  vector<double> cvalue(3*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[3*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<3>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<3>::Tensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {9*n};
  vector<double> cvalue(9*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[9*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<3>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<3>::SymTensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {6*n};
  vector<double> cvalue(6*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[6*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a vector<Dim<3>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const vector<Dim<3>::ThirdRankTensor>& value, const string pathName) {
  const string varname = this->setDir(pathName);
  const unsigned n = value.size();
  int dims[1] = {27*n};
  vector<double> cvalue(27*n);
  for (unsigned i = 0; i != value.size(); ++i) copy(value[i].begin(), value[i].end(), &cvalue[27*i]);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*cvalue.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a vector<int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<int>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  value = vector<int>(n);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value.front()) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a vector<double> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<double>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  value = vector<double>(n);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value.front()) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a vector<string> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<string>& value, const string pathName) const {
  vector<int> dim_stuff;
  string stuff;
  this->read(dim_stuff, pathName + "/dim_stuff");
  this->read(stuff, pathName + "/stuff");
  const unsigned n = dim_stuff.size();
  value = vector<string>(n);
  unsigned i = 0;
  for (unsigned k = 0; k != n; ++k) {
    value[k] = stuff.substr(i, dim_stuff[k]);
    i += dim_stuff[k];
  }

  // const string varname = this->setDir(pathName);
  // const int n = DBGetVarLength(mFilePtr, varname.c_str());
  // VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  // char** cvalue = (char**) DBGetVar(mFilePtr, varname.c_str());
  // VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  // value = vector<string>(n);
  // for (unsigned i = 0; i != n; ++i) {
  //   value[i] = string(cvalue[i]);
  // }
  // for (unsigned i = 0; i != n; ++i) {
  //   delete cvalue[i];
  // }
  // delete[] cvalue;
}

//------------------------------------------------------------------------------
// Read a vector<Dim<1>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<1>::Vector>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<1>::Vector>(n);
  for (unsigned i = 0; i != n; ++i) value[i].x(cvalue[i]);
}

//------------------------------------------------------------------------------
// Read a vector<Dim<1>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<1>::Tensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<1>::Tensor>(n);
  for (unsigned i = 0; i != n; ++i) value[i].xx(cvalue[i]);
}

//------------------------------------------------------------------------------
// Read a vector<Dim<1>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<1>::SymTensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<1>::SymTensor>(n);
  for (unsigned i = 0; i != n; ++i) value[i].xx(cvalue[i]);
}

//------------------------------------------------------------------------------
// Read a vector<Dim<1>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<1>::ThirdRankTensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<1>::ThirdRankTensor>(n);
  for (unsigned i = 0; i != n; ++i) value[i](0,0,0) = cvalue[i];
}

//------------------------------------------------------------------------------
// Read a vector<Dim<2>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<2>::Vector>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 2 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<2>::Vector>(n/2);
  for (unsigned i = 0; i != n/2; ++i) copy(&cvalue[2*i], &cvalue[2*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Read a vector<Dim<2>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<2>::Tensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 4 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<2>::Tensor>(n/4);
  for (unsigned i = 0; i != n/4; ++i) copy(&cvalue[4*i], &cvalue[4*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Read a vector<Dim<2>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<2>::SymTensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 3 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<2>::SymTensor>(n/3);
  for (unsigned i = 0; i != n/3; ++i) copy(&cvalue[3*i], &cvalue[3*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Read a vector<Dim<2>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<2>::ThirdRankTensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 8 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<2>::ThirdRankTensor>(n/8);
  for (unsigned i = 0; i != n/8; ++i) copy(&cvalue[8*i], &cvalue[8*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Read a vector<Dim<3>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<3>::Vector>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 3 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<3>::Vector>(n/3);
  for (unsigned i = 0; i != n/3; ++i) copy(&cvalue[3*i], &cvalue[3*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Read a vector<Dim<3>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<3>::Tensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 9 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<3>::Tensor>(n/9);
  for (unsigned i = 0; i != n/9; ++i) copy(&cvalue[9*i], &cvalue[9*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Read a vector<Dim<3>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<3>::SymTensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 6 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<3>::SymTensor>(n/6);
  for (unsigned i = 0; i != n/6; ++i) copy(&cvalue[6*i], &cvalue[6*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Read a vector<Dim<3>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<Dim<3>::ThirdRankTensor>& value, const string pathName) const {
  const string varname = this->setDir(pathName);
  const int n = DBGetVarLength(mFilePtr, varname.c_str());
  VERIFY2(n >= 0, "SileFileIO ERROR: unable to get size of " << pathName);
  CHECK(n % 27 == 0);
  double* cvalue = (double*) DBGetVar(mFilePtr, varname.c_str());
  VERIFY2(cvalue != NULL, "SiloFileIO Error: unable to read " << pathName);
  value = vector<Dim<3>::ThirdRankTensor>(n/27);
  for (unsigned i = 0; i != n/27; ++i) copy(&cvalue[27*i], &cvalue[27*(i+1)], value[i].begin());
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<double> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<1>::Vector> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<1>::Tensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<1>::SymTensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<1>::ThirdRankTensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, int>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<int> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<double> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<2>::Vector> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<2>::Tensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<2>::SymTensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<2>::ThirdRankTensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, int>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<int> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<double> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<3>::Vector> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<3>::Tensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<3>::SymTensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<Dim<3>::ThirdRankTensor> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, int>& value, const string pathName) {
  this->write(value.name(), pathName + "/name");
  std::vector<int> values(value.begin(), value.begin() + value.numInternalElements());
  this->write(values, pathName + "/values");
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) const {
  string fieldname;
  std::vector<double> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<1>::Vector> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<1>::Tensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<1>::SymTensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<1>::ThirdRankTensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, int>& value, const string pathName) const {
  string fieldname;
  std::vector<int> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) const {
  string fieldname;
  std::vector<double> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<2>::Vector> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<2>::Tensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<2>::SymTensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<2>::ThirdRankTensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, int>& value, const string pathName) const {
  string fieldname;
  std::vector<int> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) const {
  string fieldname;
  std::vector<double> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<3>::Vector> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<3>::Tensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<3>::SymTensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) const {
  string fieldname;
  std::vector<Dim<3>::ThirdRankTensor> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, int>& value, const string pathName) const {
  string fieldname;
  std::vector<int> values;
  this->read(fieldname, pathName + "/name");
  this->read(values, pathName + "/values");
  value.name(fieldname);
  CHECK(value.numInternalElements() == values.size());
  copy(values.begin(), values.end(), value.begin());
}

//------------------------------------------------------------------------------
// Set the directory in the current silo file to the path portion of a pathName,
// and return the variable name section.
//------------------------------------------------------------------------------
// Non-const version creates directory.
string 
SiloFileIO::setDir(const string& ipathName) {
  REQUIRE(mFileOpen and mFilePtr != 0);

  string pathName = SILO_mangle(ipathName);

  // We always start from the top.
  VERIFY2(DBSetDir(mFilePtr, "/") == 0,
          "SiloFileIO ERROR: unable to change to path /");

  // If there is no absolute path specified, we just go with the cwd.
  const size_t i = pathName.find_last_of("/");
  if (i >= pathName.size()) {
    return pathName;
  }

  // Otherwise set the path and return just the variable name.
  const string dirName = pathName.substr(0, i);
  const string varName = pathName.substr(i+1, pathName.size());

  // We have to make the directories one at a time, so progressively split
  // the path.
  if (dirName.size() > 0) {
    size_t pos0 = 0, pos1 = min(dirName.size(), dirName.find_first_of("/"));
    bool done = false;
    while (not done) {
      done = (pos1 == dirName.size());
      const string partialdir = dirName.substr(pos0, pos1-pos0);

      // Check if this directory already exists or not.
      DBtoc* toc = DBGetToc(mFilePtr);
      bool exists = false;
      size_t j = 0;
      while (not exists and j != toc->ndir) {
        exists = (partialdir == string(toc->dir_names[j]));
        ++j;
      }
      if (not exists) {
        VERIFY2(DBMkDir(mFilePtr, partialdir.c_str()) == 0,
                "SiloFileIO ERROR: unable to create path " << partialdir);
      }
      VERIFY2(DBSetDir(mFilePtr, partialdir.c_str()) == 0,
              "SiloFileIO ERROR: unable to change to path " << partialdir);

      if (not done) {
        pos0 = pos1 + 1;
        const string tail = dirName.substr(pos0, dirName.size());
        pos1 = tail.find_first_of("/");
        if (pos1 >= tail.size()) {
          pos1 = dirName.size();
        } else {
          pos1 += pos0;
        }
      }
    }
  }

  return varName;
}

// const version only changes directory.
string 
SiloFileIO::setDir(const string& ipathName) const {
  REQUIRE(mFileOpen and mFilePtr != 0);

  string pathName = SILO_mangle(ipathName);

  // We always start from the top.
  VERIFY2(DBSetDir(mFilePtr, "/") == 0,
          "SiloFileIO ERROR: unable to change to path /");

  // If there is no absolute path specified, we just go with the cwd.
  const size_t i = pathName.find_last_of("/");
  if (i >= pathName.size()) {
    return pathName;
  }

  // Otherwise set the path and return just the variable name.
  const string dirName = pathName.substr(0,i);
  const string varName = pathName.substr(i+1);
  VERIFY2(DBSetDir(mFilePtr, dirName.c_str()) == 0,
          "SiloFileIO ERROR: unable to change to path " << dirName);
  return varName;
}

}
}
