#include <iostream>
#include <vector>

#include "ExtendFileIO.hh"
#include "cdebug.hh"

namespace Spheral {
namespace FileIOSpace {

using namespace std;
using namespace boost::python;

//------------------------------------------------------------------------------
// Define helper methods to convert lists to STL vectors and call the 
// appropriate FileIO method.
// Note this assumes that the list is of a homogeneous type DataType.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
writeGenericList(FileIO& self,
                 boost::python::list& valueList,
                 const std::string pathName) {

  const std::size_t n = extract<std::size_t>(valueList.attr("__len__")());
  cdebug << "Got length of " << n << endl;
  std::vector<DataType> valueVec;
  valueVec.reserve(n);

  for (std::size_t i = 0; i < n; ++i) {
    DataType val = boost::python::extract<DataType>(valueList[i]);
//     DataType val = boost::python::extract<DataType>(valueList.attr("__getitem__")(i));
    valueVec.push_back(val);
  }

  self.write(valueVec, pathName);
}

void
writeList(FileIO& self,
          boost::python::list& valueList,
          const std::string pathName) {

  // Try to detect the element type of the list, and call the appropriate
  // instantiation of the above method.  If the thing is zero sized, though,
  // it doesn't matter which version we call.
  const std::size_t n = extract<std::size_t>(valueList.attr("__len__")());
  if (n > 0) {
    // Get the first element in the sequence.
    boost::python::object e = valueList[0];

    // Build the trial extractors.
    boost::python::extract<int> xint(e);
    boost::python::extract<float> xfloat(e);
    boost::python::extract<Spheral::Dim<1>::Vector> xvec1(e);
    boost::python::extract<Spheral::Dim<2>::Vector> xvec2(e);
    boost::python::extract<Spheral::Dim<3>::Vector> xvec3(e);
    boost::python::extract<Spheral::Dim<1>::Tensor> xtensor1(e);
    boost::python::extract<Spheral::Dim<2>::Tensor> xtensor2(e);
    boost::python::extract<Spheral::Dim<3>::Tensor> xtensor3(e);
    boost::python::extract<Spheral::Dim<1>::SymTensor> xsymtensor1(e);
    boost::python::extract<Spheral::Dim<2>::SymTensor> xsymtensor2(e);
    boost::python::extract<Spheral::Dim<3>::SymTensor> xsymtensor3(e);

    if (xint.check()) {
      cdebug << "Try list of int...";
      writeGenericList<int>(self, valueList, pathName);

    } else if (xfloat.check()) {
      cdebug << "Try list of float...";
      writeGenericList<float>(self, valueList, pathName);

    } else if (xvec1.check()) {
      cdebug << "Try list of Vector1d...";
      writeGenericList<Spheral::Dim<1>::Vector>(self, valueList, pathName);

    } else if (xvec2.check()) {
      cdebug << "Try list of Vector2d...";
      writeGenericList<Spheral::Dim<2>::Vector>(self, valueList, pathName);

    } else if (xvec3.check()) {
      cdebug << "Try list of Vector3d...";
      writeGenericList<Spheral::Dim<3>::Vector>(self, valueList, pathName);

    } else if (xtensor1.check()) {
      cdebug << "Try list of Tensor1d...";
      writeGenericList<Spheral::Dim<1>::Tensor>(self, valueList, pathName);

    } else if (xtensor2.check()) {
      cdebug << "Try list of Tensor2d...";
      writeGenericList<Spheral::Dim<2>::Tensor>(self, valueList, pathName);

    } else if (xtensor3.check()) {
      cdebug << "Try list of Tensor3d...";
      writeGenericList<Spheral::Dim<3>::Tensor>(self, valueList, pathName);

    } else if (xsymtensor1.check()) {
      cdebug << "Try list of SymTensor1d...";
      writeGenericList<Spheral::Dim<1>::SymTensor>(self, valueList, pathName);

    } else if (xsymtensor2.check()) {
      cdebug << "Try list of SymTensor2d...";
      writeGenericList<Spheral::Dim<2>::SymTensor>(self, valueList, pathName);

    } else if (xsymtensor3.check()) {
      cdebug << "Try list of SymTensor3d...";
      writeGenericList<Spheral::Dim<3>::SymTensor>(self, valueList, pathName);

    } else {
      cdebug << "ERROR: List of unknown datatype." << endl;
    }

  } else {
    // This list is empty, so just write it out as an int list and be done with it.
    writeGenericList<int>(self, valueList, pathName);
  }

}

//------------------------------------------------------------------------------
// Return a list containing elements from a std::vector write.
//------------------------------------------------------------------------------
template<typename DataType>
inline
boost::python::list
readGenericList(FileIO& self,
                const std::string pathName) {

  // Construct an empty std::vector to read the result in.
  vector<DataType> resultVector;
  self.read(resultVector, pathName);

  // Construct a list from the vector, and return that result.
  boost::python::list result(resultVector);

  return result;
}

// list of int
boost::python::list
readIntList(FileIO& self,
            const std::string pathName) {
  return readGenericList<int>(self, pathName);
}

// list of float
boost::python::list
readFloatList(FileIO& self,
              const std::string pathName) {
  return readGenericList<float>(self, pathName);
}

// list of Spheral::Dim<1>::Vector
boost::python::list
readVector1dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<1>::Vector>(self, pathName);
}

// list of Spheral::Dim<2>::Vector
boost::python::list
readVector2dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<2>::Vector>(self, pathName);
}

// list of Spheral::Dim<3>::Vector
boost::python::list
readVector3dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<3>::Vector>(self, pathName);
}

// list of Spheral::Dim<1>::Tensor
boost::python::list
readTensor1dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<1>::Tensor>(self, pathName);
}

// list of Spheral::Dim<2>::Tensor
boost::python::list
readTensor2dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<2>::Tensor>(self, pathName);
}

// list of Spheral::Dim<3>::Tensor
boost::python::list
readTensor3dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<3>::Tensor>(self, pathName);
}

// list of Spheral::Dim<1>::SymTensor
boost::python::list
readSymTensor1dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<1>::SymTensor>(self, pathName);
}

// list of Spheral::Dim<2>::SymTensor
boost::python::list
readSymTensor2dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<2>::SymTensor>(self, pathName);
}

// list of Spheral::Dim<3>::SymTensor
boost::python::list
readSymTensor3dList(FileIO& self,
            const std::string pathName) {
  return readGenericList<Spheral::Dim<3>::SymTensor>(self, pathName);
}

//------------------------------------------------------------------------------
// Force the reading of a string.
//------------------------------------------------------------------------------
std::string
readString(FileIO& self,
           const std::string pathName) {
  std::string result;
  self.read(result, pathName);
  return result;
}

}
}
