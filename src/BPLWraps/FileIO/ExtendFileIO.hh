#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Field/Field.hh"

#ifndef __GCCXML__
#include "boost/python.hpp"
#else
#include "../fakeboost.hh"
#endif

namespace Spheral {
namespace FileIOSpace {

using namespace std;
using namespace boost::python;

void
writeList(FileIO& self,
          boost::python::list& valueList,
          const std::string pathName);

// list of int
boost::python::list
readIntList(FileIO& self,
            const std::string pathName);

// list of float
boost::python::list
readFloatList(FileIO& self,
              const std::string pathName);

// list of Spheral::Dim<1>::Vector
boost::python::list
readVector1dList(FileIO& self,
                 const std::string pathName);

// list of Spheral::Dim<2>::Vector
boost::python::list
readVector2dList(FileIO& self,
                 const std::string pathName);

// list of Spheral::Dim<3>::Vector
boost::python::list
readVector3dList(FileIO& self,
                 const std::string pathName);

// list of Spheral::Dim<1>::Tensor
boost::python::list
readTensor1dList(FileIO& self,
                 const std::string pathName);

// list of Spheral::Dim<2>::Tensor
boost::python::list
readTensor2dList(FileIO& self,
                 const std::string pathName);

// list of Spheral::Dim<3>::Tensor
boost::python::list
readTensor3dList(FileIO& self,
                 const std::string pathName);

// list of Spheral::Dim<1>::SymTensor
boost::python::list
readSymTensor1dList(FileIO& self,
                    const std::string pathName);

// list of Spheral::Dim<2>::SymTensor
boost::python::list
readSymTensor2dList(FileIO& self,
                    const std::string pathName);

// list of Spheral::Dim<3>::SymTensor
boost::python::list
readSymTensor3dList(FileIO& self,
                    const std::string pathName);

//------------------------------------------------------------------------------
// Force the reading of a string.
//------------------------------------------------------------------------------
std::string
readString(FileIO& self,
           const std::string pathName);

}
}
