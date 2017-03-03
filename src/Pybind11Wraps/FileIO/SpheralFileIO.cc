// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "FileIO/FileIO.hh"
#include "FileIO/FlatFileIO.hh"
#include "FileIO/SiloFileIO.hh"
#include "FileIO/PyFileIO.hh"
#include "FileIO/vectorstringUtilities.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::FileIOSpace;
using namespace Spheral::FieldSpace;

namespace Spheral {
namespace FileIOSpace {

//------------------------------------------------------------------------------
// PyAbstractFileIO
//------------------------------------------------------------------------------
template<class FileIOBase>
class PyAbstractFileIO: public FileIOBase {
public:
  using FileIOBase::FileIOBase;  // inherit constructors

  virtual void open(const std::string fileName, AccessType access) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, open, fileName, access); }
  virtual void close() override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, close); }

  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const unsigned value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const int value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const bool value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const double value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::string value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void write(const Dim<1>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<1>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<1>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void write(const Dim<2>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<2>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<2>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void write(const Dim<3>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<3>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<3>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void read(unsigned& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(int& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(bool& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(double& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::string& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  virtual void read(Dim<1>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<1>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<1>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  virtual void read(Dim<2>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<2>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<2>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  virtual void read(Dim<3>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<3>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<3>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  // We also require that FileIO objects write vectors of the primitive types.
  virtual void write(const std::vector<int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<double>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<std::string>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void write(const std::vector<Dim<1>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<1>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<1>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void write(const std::vector<Dim<2>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<2>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<2>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void write(const std::vector<Dim<3>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<3>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<3>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void read(std::vector<int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<double>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<std::string>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  virtual void read(std::vector<Dim<1>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<1>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<1>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  virtual void read(std::vector<Dim<2>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<2>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<2>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  virtual void read(std::vector<Dim<3>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<3>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<3>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }

  // Require that all FileIO objects provide methods to read and write
  // Fields of specific DataTypes.
#ifdef SPHERAL1D
  virtual void write(const Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void read(Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
#endif

#ifdef SPHERAL2D
  virtual void write(const Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void read(Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
#endif

#ifdef SPHERAL3D
  virtual void write(const Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write, value, pathName); }

  virtual void read(Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read, value, pathName); }
#endif

  // These methods are useful for the primitive types that are problematic
  // to return by reference from python.
  virtual void write_unsigned_int(const unsigned value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write_unsigned_int, value, pathName); }
  virtual void write_int(const int value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write_int, value, pathName); }
  virtual void write_bool(const bool value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write_bool, value, pathName); }
  virtual void write_double(const double value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write_double, value, pathName); }
  virtual void write_string(const std::string value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write_string, value, pathName); }

  virtual unsigned    read_unsigned_int(const std::string pathName) const override { PYBIND11_OVERLOAD(unsigned, FileIOBase, read_unsigned_int, pathName); }
  virtual int         read_int(const std::string pathName) const override { PYBIND11_OVERLOAD(int, FileIOBase, read_int, pathName); }
  virtual bool        read_bool(const std::string pathName) const override { PYBIND11_OVERLOAD(bool, FileIOBase, read_bool, pathName); }
  virtual double      read_double(const std::string pathName) const override { PYBIND11_OVERLOAD(double, FileIOBase, read_double, pathName); }
  virtual std::string read_string(const std::string pathName) const override { PYBIND11_OVERLOAD(std::string, FileIOBase, read_string, pathName); }
};

}
}

namespace {  // anonymous

//------------------------------------------------------------------------------
// Common virtual methods of FileIO objects.
//------------------------------------------------------------------------------
template<typename Obj, typename PB11Obj>
void virtualFileIOBindings(py::module& m, PB11Obj& obj) {

  obj
    .def("write", (void (Obj::*)(const unsigned, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const int, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const bool, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const double, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::string, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    ;
}

}

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralFileIO) {
  py::module m("SpheralFileIO", "Spheral FileIO module.");

  //............................................................................
  // AccessType
  py::enum_<Spheral::FileIOSpace::AccessType>(m, "AccessType")
    .value("Undefined", Spheral::FileIOSpace::AccessType::Undefined)
    .value("Create", Spheral::FileIOSpace::AccessType::Create)
    .value("Read", Spheral::FileIOSpace::AccessType::Read)
    .value("Write", Spheral::FileIOSpace::AccessType::Write)
    .value("ReadWrite", Spheral::FileIOSpace::AccessType::ReadWrite)
    .export_values();

  //............................................................................
  // FlatFileFormat
  py::enum_<Spheral::FileIOSpace::FlatFileFormat>(m, "FlatFileFormat")
    .value("ascii", Spheral::FileIOSpace::FlatFileFormat::ascii)
    .value("binary", Spheral::FileIOSpace::FlatFileFormat::binary)
    .export_values();

  //............................................................................
  // FileIO
  py::class_<Spheral::FileIOSpace::FileIO, PyAbstractFileIO<Spheral::FileIOSpace::FileIO>> fileioPB11(m, "FileIO");
  virtualFileIOBindings<Spheral::FileIOSpace::FileIO>(m, fileioPB11);

  fileioPB11

    // Constructors
    .def(py::init<>())
    .def(py::init<std::string, Spheral::FileIOSpace::AccessType>(), "filename"_a, "access"_a)

    ;

  return m.ptr();
}
