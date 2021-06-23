//---------------------------------Spheral++----------------------------------//
// SiloFileIO -- Provide the interface to silo file objects.
//
// Created by JMO, Sat Feb  7 23:06:03 PST 2015
//----------------------------------------------------------------------------//
#include "SiloFileIO.hh"
#include "Field/Field.hh"

#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/replace.hpp"

#include <algorithm>
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

namespace {
  
//------------------------------------------------------------------------------
// Mangle a path name for silo consumption.
//------------------------------------------------------------------------------
string SILO_mangle(const string& x) {
  string result = x;
  boost::replace_all(result, " ", "_");
  boost::replace_all(result, ".", "_p_");
  boost::replace_all(result, ",", "_c_");
  boost::replace_all(result, "|", "_P_");
  boost::replace_all(result, "-", "_d_");
  return result;
}

//------------------------------------------------------------------------------
// Set the directory in the current silo file to the path portion of a pathName,
// and return the variable name section.
//------------------------------------------------------------------------------
string setdir(DBfile* filePtr, const string& ipathName) {

  // Mangle the path into something acceptable to silo.
  string pathName = SILO_mangle(ipathName);

  // We always start from the top.
  VERIFY2(DBSetDir(filePtr, "/") == 0,
          "SiloFileIO ERROR: unable to change to path /");

  // Split up the path into directories and the var name.
  vector<string> components;
  boost::split(components, pathName, boost::is_any_of("/"));
  CHECK(components.size() > 0);

  // Set the path and return just the variable name.
  for (unsigned i = 0; i != components.size() - 1; ++i) {
    const string dirName = components[i];
    if (dirName.size() > 0) {
      // Check if this directory already exists or not.
      DBtoc* toc = DBGetToc(filePtr);
      bool exists = false;
      int j = 0;
      while (not exists and j != toc->ndir) {
        exists = (dirName == string(toc->dir_names[j]));
        ++j;
      }
      if (not exists) {
        // std::cerr << " --> " << pathName << std::endl
        //           << " --> " << dirName.size() << std::endl
        //           << " --> " << dirName << std::endl;
        VERIFY2(DBMkDir(filePtr, dirName.c_str()) == 0,
                "SiloFileIO ERROR: unable to create path " << dirName << " of " << pathName);
      }
      VERIFY2(DBSetDir(filePtr, dirName.c_str()) == 0,
              "SiloFileIO ERROR: unable to change to path " << dirName << " of " << pathName);
    }
  }

  return components.back();
}


//------------------------------------------------------------------------------
// Common methods for reading/writing a lowly int to a silo file.
//------------------------------------------------------------------------------
void writeInt(DBfile* filePtr, const int& value, const string pathName) {
  const string varname = setdir(filePtr, pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(filePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write int variable " << pathName);
}

void readInt(DBfile* filePtr, int& value, const string pathName) {
  // const string varname = setdir(mFilePtr, pathName);
  // CHECK2(DBReadVar(filePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
  const string varname = setdir(filePtr, pathName);
  value = *static_cast<int*>(DBGetVar(filePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Common methods for reading/writing a std::string to a silo file.
//------------------------------------------------------------------------------
void writeString(DBfile* filePtr, const string& value, const string pathName) {
  const int size = value.size();
  char cvalue[size];
  std::copy(value.begin(), value.end(), cvalue);
  //cvalue[size] = '\0';
  int dims[1] = {size};
  writeInt(filePtr, dims[0], pathName + "/size");
  if (dims[0] > 0) {
    const string varname = setdir(filePtr, pathName + "/value");
    VERIFY2(DBWrite(filePtr, varname.c_str(), cvalue, dims, 1, DB_CHAR) == 0,
            "SiloFileIO ERROR: unable to write string variable " << pathName);
  }
}

void readString(DBfile* filePtr, string& value, const string pathName) {
  int valsize;
  readInt(filePtr, valsize, pathName + "/size");
  if (valsize == 0) {
    value = "";
  } else {
    char cvalue[valsize + 1];  // Do we need to allow space for trailing null?
    const string varname = setdir(filePtr, pathName + "/value");
    VERIFY2(DBReadVar(filePtr, varname.c_str(), cvalue) == 0,
            "SiloFileIO ERROR: failed to read string variable " << pathName);
    value = string(&cvalue[0], &cvalue[valsize]);
  }
}

//------------------------------------------------------------------------------
// Generic methods to read/write a Field
// These methods work for types Vector, Tensor, ThirdRankTensor, ...
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
struct FieldIO {

  //...........................................................................
  static void write(DBfile* filePtr,
                    const Field<Dimension, Value>& value, 
                    const string pathName) {
    const int n = value.numInternalElements();
    const int ne = Value::numElements;
    int dims[1] = {n*ne};
    writeInt(filePtr, dims[0], pathName + "/values/size");
    writeString(filePtr, value.name(), pathName + "/name");
    if (n > 0) {
      std::vector<double> buf(n*ne);
      for (auto i = 0; i < n; ++i) std::copy(value(i).begin(), value(i).end(), &buf[i*ne]);
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBWrite(filePtr, varname.c_str(), static_cast<void*>(&buf.front()), dims, 1, DB_DOUBLE) == 0,
              "SiloFileIO ERROR: unable to write Field values " << pathName);
    }
  }

  //...........................................................................
  static void read(DBfile* filePtr,
                   Field<Dimension, Value>& value, 
                   const string pathName) {
    const int n = value.numInternalElements();
    const int ne = Value::numElements;
    int n1;
    string fieldname;
    readInt(filePtr, n1, pathName + "/values/size");
    readString(filePtr, fieldname, pathName + "/name");
    VERIFY2(n*ne == n1, "SiloFileIO ERROR: bad Field size " << n << " != " << n1);
    value.name(fieldname);
    if (n > 0) {
      std::vector<double> buf(n*ne);
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBReadVar(filePtr, varname.c_str(), static_cast<void*>(&buf.front())) == 0,
              "SiloFileIO ERROR: unable to read Field values " << pathName);
      for (auto i = 0; i < n; ++i) std::copy(&buf[i*ne], &buf[i*ne] + ne, value(i).begin());
    }
  }
};

//------------------------------------------------------------------------------
// Specialize to read/write a Field<int>
//------------------------------------------------------------------------------
template<typename Dimension>
struct FieldIO<Dimension, int> {

  //...........................................................................
  static void write(DBfile* filePtr,
                    const Field<Dimension, int>& value, 
                    const string pathName) {
    const int n = value.numInternalElements();
    writeInt(filePtr, n, pathName + "/values/size");
    writeString(filePtr, value.name(), pathName + "/name");
    if (n > 0) {
      int dims[1] = {n};
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBWrite(filePtr, varname.c_str(), static_cast<void*>(const_cast<int*>(&value[0])), dims, 1, DB_INT) == 0,
              "SiloFileIO ERROR: unable to write Field values " << pathName);
    }
  }

  //...........................................................................
  static void read(DBfile* filePtr,
                   Field<Dimension, int>& value, 
                   const string pathName) {
    const int n = value.numInternalElements();
    int n1;
    string fieldname;
    readInt(filePtr, n1, pathName + "/values/size");
    readString(filePtr, fieldname, pathName + "/name");
    VERIFY2(n == n1, "SiloFileIO ERROR: bad Field size " << n << " != " << n1);
    value.name(fieldname);
    if (n > 0) {
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBReadVar(filePtr, varname.c_str(), static_cast<void*>(&value[0])) == 0,
              "SiloFileIO ERROR: unable to read Field values " << pathName);
    }
  }
};


//------------------------------------------------------------------------------
// Specialize to read/write a Field<unsigned>
//------------------------------------------------------------------------------
template<typename Dimension>
struct FieldIO<Dimension, unsigned> {

  //...........................................................................
  static void write(DBfile* filePtr,
                    const Field<Dimension, unsigned>& value, 
                    const string pathName) {
    const int n = value.numInternalElements();
    writeInt(filePtr, n, pathName + "/values/size");
    writeString(filePtr, value.name(), pathName + "/name");
    if (n > 0) {
      int dims[1] = {n};
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBWrite(filePtr, varname.c_str(), static_cast<void*>(const_cast<unsigned*>(&value[0])), dims, 1, DB_INT) == 0,
              "SiloFileIO ERROR: unable to write Field values " << pathName);
    }
  }

  //...........................................................................
  static void read(DBfile* filePtr,
                   Field<Dimension, unsigned>& value, 
                   const string pathName) {
    const int n = value.numInternalElements();
    int n1;
    string fieldname;
    readInt(filePtr, n1, pathName + "/values/size");
    readString(filePtr, fieldname, pathName + "/name");
    VERIFY2(n == n1, "SiloFileIO ERROR: bad Field size " << n << " != " << n1);
    value.name(fieldname);
    if (n > 0) {
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBReadVar(filePtr, varname.c_str(), static_cast<void*>(&value[0])) == 0,
              "SiloFileIO ERROR: unable to read Field values " << pathName);
    }
  }
};


//------------------------------------------------------------------------------
// Specialize to read/write a Field<double>
//------------------------------------------------------------------------------
template<typename Dimension>
struct FieldIO<Dimension, typename Dimension::Scalar> {

  //...........................................................................
  static void write(DBfile* filePtr,
                    const Field<Dimension, double>& value, 
                    const string pathName) {
    const int n = value.numInternalElements();
    writeInt(filePtr, n, pathName + "/values/size");
    writeString(filePtr, value.name(), pathName + "/name");
    if (n > 0) {
      int dims[1] = {n};
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBWrite(filePtr, varname.c_str(), static_cast<void*>(const_cast<double*>(&value[0])), dims, 1, DB_DOUBLE) == 0,
              "SiloFileIO ERROR: unable to write Field values " << pathName);
    }
  }

  //...........................................................................
  static void read(DBfile* filePtr,
                   Field<Dimension, double>& value, 
                   const string pathName) {
    const int n = value.numInternalElements();
    int n1;
    string fieldname;
    readInt(filePtr, n1, pathName + "/values/size");
    readString(filePtr, fieldname, pathName + "/name");
    VERIFY2(n == n1, "SiloFileIO ERROR: bad Field size " << n << " != " << n1);
    value.name(fieldname);
    if (n > 0) {
      const string varname = setdir(filePtr, pathName + "/values/value");
      VERIFY2(DBReadVar(filePtr, varname.c_str(), static_cast<void*>(&value[0])) == 0,
              "SiloFileIO ERROR: unable to read Field values " << pathName);
    }
  }
};

}   // anonymous namespace

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

  if (access == AccessType::Read) {
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
// Check if the specified path is in the file.
//------------------------------------------------------------------------------
bool
SiloFileIO::pathExists(const std::string pathName) const {
  return DBInqVarExists(mFilePtr, pathName.c_str()) != 0;
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const unsigned& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a size_t to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const size_t& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const int& value, const string pathName) {
  writeInt(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const bool& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {1};
  int ivalue = value ? 1 : 0;
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &ivalue, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const double& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {1};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const string& value, const string pathName) {
  writeString(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a vector<int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const std::vector<int>& value, const string pathName) {
  const int size = value.size();
  writeInt(mFilePtr, size, pathName + "/size");
  if (size > 0) {
    const string varname = setdir(mFilePtr, pathName + "/value");
    int dims[1] = {size};
    VERIFY2(DBWrite(mFilePtr, varname.c_str(), static_cast<void*>(const_cast<int*>(&value[0])), dims, 1, DB_INT) == 0,
            "SiloFileIO ERROR: unable to write std::vector " << pathName);
  }
}

//------------------------------------------------------------------------------
// Write a vector<double> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const std::vector<double>& value, const string pathName) {
  const int size = value.size();
  writeInt(mFilePtr, size, pathName + "/size");
  if (size > 0) {
    const string varname = setdir(mFilePtr, pathName + "/value");
    int dims[1] = {size};
    VERIFY2(DBWrite(mFilePtr, varname.c_str(), static_cast<void*>(const_cast<double*>(&value[0])), dims, 1, DB_DOUBLE) == 0,
            "SiloFileIO ERROR: unable to write std::vector " << pathName);
  }
}

//------------------------------------------------------------------------------
// Write a vector<string> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const std::vector<string>& value, const string pathName) {
  const unsigned n = value.size();
  vector<int> dim_stuff(n);
  string stuff;
  for (unsigned i = 0; i != n; ++i) {
    dim_stuff[i] = value[i].size();
    stuff += value[i];
  }
  this->write(dim_stuff, pathName + "/dim_stuff");
  this->write(stuff, pathName + "/stuff");
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Vector& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<1>::Vector::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Tensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<1>::Tensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::SymTensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<1>::SymTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::ThirdRankTensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<1>::ThirdRankTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Vector& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<2>::Vector::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Tensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<2>::Tensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::SymTensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<2>::SymTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::ThirdRankTensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<2>::ThirdRankTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Vector& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<3>::Vector::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Tensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<3>::Tensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::SymTensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<3>::SymTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::ThirdRankTensor& value, const string pathName) {
  const string varname = setdir(mFilePtr, pathName);
  int dims[1] = {int(Dim<3>::ThirdRankTensor::numElements)};
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << pathName);
}

//------------------------------------------------------------------------------
// Read an unsigned from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(unsigned& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
  value = *static_cast<unsigned*>(DBGetVar(mFilePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Read a size_t from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(size_t& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
  value = *static_cast<size_t*>(DBGetVar(mFilePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Read an int to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(int& value, const string pathName) const {
  readInt(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(bool& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  const int ivalue = *static_cast<int*>(DBGetVar(mFilePtr, varname.c_str()));
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &ivalue) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
  value = (ivalue == 1 ? true : false);
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(double& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  value = *static_cast<double*>(DBGetVar(mFilePtr, varname.c_str()));
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(string& value, const string pathName) const {
  readString(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a vector<int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(std::vector<int>& value, const string pathName) const {
  int size;
  readInt(mFilePtr, size, pathName + "/size");
  value.resize(size);
  if (size > 0) {
    const string varname = setdir(mFilePtr, pathName + "/value");
    VERIFY2(DBReadVar(mFilePtr, varname.c_str(), static_cast<void*>(&value[0])) == 0,
            "SiloFileIO ERROR: unable to read std::vector " << pathName);
  }
}

//------------------------------------------------------------------------------
// Read a vector<double> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(std::vector<double>& value, const string pathName) const {
  int size;
  readInt(mFilePtr, size, pathName + "/size");
  value.resize(size);
  if (size > 0) {
    const string varname = setdir(mFilePtr, pathName + "/value");
    VERIFY2(DBReadVar(mFilePtr, varname.c_str(), static_cast<void*>(&value[0])) == 0,
            "SiloFileIO ERROR: unable to read std::vector " << pathName);
  }
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
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Vector& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Tensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::SymTensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::ThirdRankTensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Vector& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Tensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::SymTensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::ThirdRankTensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Vector& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Tensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::SymTensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::ThirdRankTensor& value, const string pathName) const {
  const string varname = setdir(mFilePtr, pathName);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << pathName);
}

#ifdef SPHERAL1D
//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) {
  FieldIO<Dim<1>, Dim<1>::Scalar>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) {
  FieldIO<Dim<1>, Dim<1>::Vector>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) {
  FieldIO<Dim<1>, Dim<1>::Tensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) {
  FieldIO<Dim<1>, Dim<1>::SymTensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, Dim<1>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) {
  FieldIO<Dim<1>, Dim<1>::ThirdRankTensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, int>& value, const string pathName) {
  FieldIO<Dim<1>, int>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<1>, unsigned> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<1>, unsigned>& value, const string pathName) {
  FieldIO<Dim<1>, unsigned>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Scalar>& value, const string pathName) const {
  FieldIO<Dim<1>, Dim<1>::Scalar>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Vector>& value, const string pathName) const {
  FieldIO<Dim<1>, Dim<1>::Vector>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::Tensor>& value, const string pathName) const {
  FieldIO<Dim<1>, Dim<1>::Tensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::SymTensor>& value, const string pathName) const {
  FieldIO<Dim<1>, Dim<1>::SymTensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, Dim<1>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const string pathName) const {
  FieldIO<Dim<1>, Dim<1>::ThirdRankTensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, int>& value, const string pathName) const {
  FieldIO<Dim<1>, int>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<1>, unsigned> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<1>, unsigned>& value, const string pathName) const {
  FieldIO<Dim<1>, unsigned>::read(mFilePtr, value, pathName);
}
#endif

#ifdef SPHERAL2D
//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) {
  FieldIO<Dim<2>, Dim<2>::Scalar>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) {
  FieldIO<Dim<2>, Dim<2>::Vector>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) {
  FieldIO<Dim<2>, Dim<2>::Tensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) {
  FieldIO<Dim<2>, Dim<2>::SymTensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, Dim<2>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) {
  FieldIO<Dim<2>, Dim<2>::ThirdRankTensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, int>& value, const string pathName) {
  FieldIO<Dim<2>, int>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<2>, unsigned> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<2>, unsigned>& value, const string pathName) {
  FieldIO<Dim<2>, unsigned>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Scalar>& value, const string pathName) const {
  FieldIO<Dim<2>, Dim<2>::Scalar>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Vector>& value, const string pathName) const {
  FieldIO<Dim<2>, Dim<2>::Vector>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::Tensor>& value, const string pathName) const {
  FieldIO<Dim<2>, Dim<2>::Tensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::SymTensor>& value, const string pathName) const {
  FieldIO<Dim<2>, Dim<2>::SymTensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, Dim<2>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const string pathName) const {
  FieldIO<Dim<2>, Dim<2>::ThirdRankTensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, int>& value, const string pathName) const {
  FieldIO<Dim<2>, int>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<2>, unsigned> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<2>, unsigned>& value, const string pathName) const {
  FieldIO<Dim<2>, unsigned>::read(mFilePtr, value, pathName);
}
#endif

#ifdef SPHERAL3D
//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Scalar> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) {
  FieldIO<Dim<3>, Dim<3>::Scalar>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Vector> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) {
  FieldIO<Dim<3>, Dim<3>::Vector>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::Tensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) {
  FieldIO<Dim<3>, Dim<3>::Tensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::SymTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) {
  FieldIO<Dim<3>, Dim<3>::SymTensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, Dim<3>::ThirdRankTensor> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) {
  FieldIO<Dim<3>, Dim<3>::ThirdRankTensor>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, int>& value, const string pathName) {
  FieldIO<Dim<3>, int>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Write a Field<Dim<3>, unsigned> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Field<Dim<3>, unsigned>& value, const string pathName) {
  FieldIO<Dim<3>, unsigned>::write(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Scalar> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Scalar>& value, const string pathName) const {
  FieldIO<Dim<3>, Dim<3>::Scalar>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Vector> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Vector>& value, const string pathName) const {
  FieldIO<Dim<3>, Dim<3>::Vector>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::Tensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::Tensor>& value, const string pathName) const {
  FieldIO<Dim<3>, Dim<3>::Tensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::SymTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::SymTensor>& value, const string pathName) const {
  FieldIO<Dim<3>, Dim<3>::SymTensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, Dim<3>::ThirdRankTensor> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const string pathName) const {
  FieldIO<Dim<3>, Dim<3>::ThirdRankTensor>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, int>& value, const string pathName) const {
  FieldIO<Dim<3>, int>::read(mFilePtr, value, pathName);
}

//------------------------------------------------------------------------------
// Read a Field<Dim<3>, unsigned> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Field<Dim<3>, unsigned>& value, const string pathName) const {
  FieldIO<Dim<3>, unsigned>::read(mFilePtr, value, pathName);
}
#endif

}
