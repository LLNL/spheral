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

using namespace Spheral;
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

//------------------------------------------------------------------------------
// PyConcreteFileIO
//------------------------------------------------------------------------------
template<class FileIOBase>
class PyConcreteFileIO: public FileIOBase {
public:
  using FileIOBase::FileIOBase;  // inherit constructors

  virtual void open(const std::string fileName, AccessType access) override { PYBIND11_OVERLOAD(void, FileIOBase, open, fileName, access); }
  virtual void close() override { PYBIND11_OVERLOAD(void, FileIOBase, close); }

  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const unsigned value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const int value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const bool value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const double value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::string value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void write(const Dim<1>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<1>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<1>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void write(const Dim<2>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<2>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<2>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void write(const Dim<3>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<3>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<3>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void read(unsigned& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(int& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(bool& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(double& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::string& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  virtual void read(Dim<1>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<1>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<1>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  virtual void read(Dim<2>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<2>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<2>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  virtual void read(Dim<3>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<3>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<3>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  // We also require that FileIO objects write vectors of the primitive types.
  virtual void write(const std::vector<int>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<double>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<std::string>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void write(const std::vector<Dim<1>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<1>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<1>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void write(const std::vector<Dim<2>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<2>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<2>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void write(const std::vector<Dim<3>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<3>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<3>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void read(std::vector<int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<double>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<std::string>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  virtual void read(std::vector<Dim<1>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<1>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<1>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  virtual void read(std::vector<Dim<2>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<2>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<2>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  virtual void read(std::vector<Dim<3>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<3>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<3>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }

  // Require that all FileIO objects provide methods to read and write
  // Fields of specific DataTypes.
#ifdef SPHERAL1D
  virtual void write(const Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<1>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void read(Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<1>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
#endif

#ifdef SPHERAL2D
  virtual void write(const Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<2>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void read(Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<2>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
#endif

#ifdef SPHERAL3D
  virtual void write(const Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }
  virtual void write(const Field<Dim<3>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD(void, FileIOBase, write, value, pathName); }

  virtual void read(Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
  virtual void read(Field<Dim<3>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD(void, FileIOBase, read, value, pathName); }
#endif

};

//------------------------------------------------------------------------------
// PyPyConcreteFileIO
//------------------------------------------------------------------------------
template<class FileIOBase>
class PyPyConcreteFileIO: public PyConcreteFileIO<FileIOBase> {
public:
  using PyConcreteFileIO<FileIOBase>::PyConcreteFileIO;  // inherit constructors

  virtual void write_Vector1d(const Dim<1>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_Vector1d, value, pathName); }
  virtual void write_Tensor1d(const Dim<1>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_Tenor1d, value, pathName); }
  virtual void write_SymTensor1d(const Dim<1>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_SymTensor1d, value, pathName); }
  virtual void write_ThirdRankTensor1d(const Dim<1>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ThirdRankTensor1d, value, pathName); }

  virtual void write_Vector2d(const Dim<2>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_Vector2d, value, pathName); }
  virtual void write_Tensor2d(const Dim<2>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_Tenor2d, value, pathName); }
  virtual void write_SymTensor2d(const Dim<2>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_SymTensor2d, value, pathName); }
  virtual void write_ThirdRankTensor2d(const Dim<2>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ThirdRankTensor2d, value, pathName); }

  virtual void write_Vector3d(const Dim<3>::Vector& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_Vector3d, value, pathName); }
  virtual void write_Tensor3d(const Dim<3>::Tensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_Tenor3d, value, pathName); }
  virtual void write_SymTensor3d(const Dim<3>::SymTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_SymTensor3d, value, pathName); }
  virtual void write_ThirdRankTensor3d(const Dim<3>::ThirdRankTensor& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ThirdRankTensor3d, value, pathName); }

  virtual void write_vector_of_int(const std::vector<int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_int, value, pathName); }
  virtual void write_vector_of_double(const std::vector<double>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_double, value, pathName); }
  virtual void write_vector_of_string(const std::vector<std::string>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_string, value, pathName); }

  virtual void write_vector_of_Vector1d(const std::vector<Dim<1>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_Vector1d, value, pathName); }
  virtual void write_vector_of_Tensor1d(const std::vector<Dim<1>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_Tensor1d, value, pathName); }
  virtual void write_vector_of_SymTensor1d(const std::vector<Dim<1>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_SymTensor1d, value, pathName); }
  virtual void write_vector_of_ThirdRankTensor1d(const std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_ThirdRankTensor1d, value, pathName); }

  virtual void write_vector_of_Vector2d(const std::vector<Dim<2>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_Vector2d, value, pathName); }
  virtual void write_vector_of_Tensor2d(const std::vector<Dim<2>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_Tensor2d, value, pathName); }
  virtual void write_vector_of_SymTensor2d(const std::vector<Dim<2>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_SymTensor2d, value, pathName); }
  virtual void write_vector_of_ThirdRankTensor2d(const std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_ThirdRankTensor2d, value, pathName); }

  virtual void write_vector_of_Vector3d(const std::vector<Dim<3>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_Vector3d, value, pathName); }
  virtual void write_vector_of_Tensor3d(const std::vector<Dim<3>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_Tensor3d, value, pathName); }
  virtual void write_vector_of_SymTensor3d(const std::vector<Dim<3>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_SymTensor3d, value, pathName); }
  virtual void write_vector_of_ThirdRankTensor3d(const std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_vector_of_ThirdRankTensor3d, value, pathName); }

#ifdef SPHERAL1D
  virtual void write_ScalarField1d(const Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ScalarField1d, value, pathName); }
  virtual void write_VectorField1d(const Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_VectorField1d, value, pathName); }
  virtual void write_TensorField1d(const Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_TensorField1d, value, pathName); }
  virtual void write_SymTensorField1d(const Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_SymTensorField1d, value, pathName); }
  virtual void write_ThirdRankTensorField1d(const Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ThirdRankTensorField1d, value, pathName); }
  virtual void write_IntField1d(const Field<Dim<1>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_IntField1d, value, pathName); }
#endif

#ifdef SPHERAL2D
  virtual void write_ScalarField2d(const Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ScalarField2d, value, pathName); }
  virtual void write_VectorField2d(const Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_VectorField2d, value, pathName); }
  virtual void write_TensorField2d(const Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_TensorField2d, value, pathName); }
  virtual void write_SymTensorField2d(const Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_SymTensorField2d, value, pathName); }
  virtual void write_ThirdRankTensorField2d(const Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ThirdRankTensorField2d, value, pathName); }
  virtual void write_IntField2d(const Field<Dim<2>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_IntField2d, value, pathName); }
#endif

#ifdef SPHERAL3D
  virtual void write_ScalarField3d(const Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ScalarField3d, value, pathName); }
  virtual void write_VectorField3d(const Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_VectorField3d, value, pathName); }
  virtual void write_TensorField3d(const Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_TensorField3d, value, pathName); }
  virtual void write_SymTensorField3d(const Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_SymTensorField3d, value, pathName); }
  virtual void write_ThirdRankTensorField3d(const Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_ThirdRankTensorField3d, value, pathName); }
  virtual void write_IntField3d(const Field<Dim<3>, int>& value, const std::string pathName) override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, write_IntField3d, value, pathName); }
#endif

  virtual void read_Vector1d(Dim<1>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_Vector1d, value, pathName); }
  virtual void read_Tensor1d(Dim<1>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_Tensor1d, value, pathName); }
  virtual void read_SymTensor1d(Dim<1>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_SymTensor1d, value, pathName); }
  virtual void read_ThirdRankTensor1d(Dim<1>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ThirdRankTensor1d, value, pathName); }

  virtual void read_Vector2d(Dim<2>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_Vector2d, value, pathName); }
  virtual void read_Tensor2d(Dim<2>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_Tensor2d, value, pathName); }
  virtual void read_SymTensor2d(Dim<2>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_SymTensor2d, value, pathName); }
  virtual void read_ThirdRankTensor2d(Dim<2>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ThirdRankTensor2d, value, pathName); }

  virtual void read_Vector3d(Dim<3>::Vector& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_Vector3d, value, pathName); }
  virtual void read_Tensor3d(Dim<3>::Tensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_Tensor3d, value, pathName); }
  virtual void read_SymTensor3d(Dim<3>::SymTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_SymTensor3d, value, pathName); }
  virtual void read_ThirdRankTensor3d(Dim<3>::ThirdRankTensor& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ThirdRankTensor3d, value, pathName); }

  virtual void read_vector_of_int(std::vector<int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_int, value, pathName); }
  virtual void read_vector_of_double(std::vector<double>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_double, value, pathName); }
  virtual void read_vector_of_string(std::vector<std::string>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_string, value, pathName); }

  virtual void read_vector_of_Vector1d(std::vector<Dim<1>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_Vector1d, value, pathName); }
  virtual void read_vector_of_Tensor1d(std::vector<Dim<1>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_Tensor1d, value, pathName); }
  virtual void read_vector_of_SymTensor1d(std::vector<Dim<1>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_SymTensor1d, value, pathName); }
  virtual void read_vector_of_ThirdRankTensor1d(std::vector<Dim<1>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_ThirdRankTensor1d, value, pathName); }

  virtual void read_vector_of_Vector2d(std::vector<Dim<2>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_Vector2d, value, pathName); }
  virtual void read_vector_of_Tensor2d(std::vector<Dim<2>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_Tensor2d, value, pathName); }
  virtual void read_vector_of_SymTensor2d(std::vector<Dim<2>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_SymTensor2d, value, pathName); }
  virtual void read_vector_of_ThirdRankTensor2d(std::vector<Dim<2>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_ThirdRankTensor2d, value, pathName); }

  virtual void read_vector_of_Vector3d(std::vector<Dim<3>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_Vector3d, value, pathName); }
  virtual void read_vector_of_Tensor3d(std::vector<Dim<3>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_Tensor3d, value, pathName); }
  virtual void read_vector_of_SymTensor3d(std::vector<Dim<3>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_SymTensor3d, value, pathName); }
  virtual void read_vector_of_ThirdRankTensor3d(std::vector<Dim<3>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_vector_of_ThirdRankTensor3d, value, pathName); }

#ifdef SPHERAL1D
  virtual void read_ScalarField1d(Field<Dim<1>, Dim<1>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ScalarField1d, value, pathName); }
  virtual void read_VectorField1d(Field<Dim<1>, Dim<1>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_VectorField1d, value, pathName); }
  virtual void read_TensorField1d(Field<Dim<1>, Dim<1>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_TensorField1d, value, pathName); }
  virtual void read_SymTensorField1d(Field<Dim<1>, Dim<1>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_SymTensorField1d, value, pathName); }
  virtual void read_ThirdRankTensorField1d(Field<Dim<1>, Dim<1>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ThirdRankTensorField1d, value, pathName); }
  virtual void read_IntField1d(Field<Dim<1>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_IntField1d, value, pathName); }
#endif

#ifdef SPHERAL2D
  virtual void read_ScalarField2d(Field<Dim<2>, Dim<2>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ScalarField2d, value, pathName); }
  virtual void read_VectorField2d(Field<Dim<2>, Dim<2>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_VectorField2d, value, pathName); }
  virtual void read_TensorField2d(Field<Dim<2>, Dim<2>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_TensorField2d, value, pathName); }
  virtual void read_SymTensorField2d(Field<Dim<2>, Dim<2>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_SymTensorField2d, value, pathName); }
  virtual void read_ThirdRankTensorField2d(Field<Dim<2>, Dim<2>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ThirdRankTensorField2d, value, pathName); }
  virtual void read_IntField2d(Field<Dim<2>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_IntField2d, value, pathName); }
#endif

#ifdef SPHERAL3D
  virtual void read_ScalarField3d(Field<Dim<3>, Dim<3>::Scalar>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ScalarField3d, value, pathName); }
  virtual void read_VectorField3d(Field<Dim<3>, Dim<3>::Vector>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_VectorField3d, value, pathName); }
  virtual void read_TensorField3d(Field<Dim<3>, Dim<3>::Tensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_TensorField3d, value, pathName); }
  virtual void read_SymTensorField3d(Field<Dim<3>, Dim<3>::SymTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_SymTensorField3d, value, pathName); }
  virtual void read_ThirdRankTensorField3d(Field<Dim<3>, Dim<3>::ThirdRankTensor>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_ThirdRankTensorField3d, value, pathName); }
  virtual void read_IntField3d(Field<Dim<3>, int>& value, const std::string pathName) const override { PYBIND11_OVERLOAD_PURE(void, FileIOBase, read_IntField3d, value, pathName); }
#endif
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

    .def("write", (void (Obj::*)(const Dim<1>::Vector&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<1>::Tensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<1>::SymTensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<1>::ThirdRankTensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    .def("write", (void (Obj::*)(const Dim<2>::Vector&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<2>::Tensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<2>::SymTensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<2>::ThirdRankTensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    .def("write", (void (Obj::*)(const Dim<3>::Vector&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<3>::Tensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<3>::SymTensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const Dim<3>::ThirdRankTensor&, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(unsigned&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(int&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(bool&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(double&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::string&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(Dim<1>::Vector&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<1>::Tensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<1>::SymTensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<1>::ThirdRankTensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(Dim<2>::Vector&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<2>::Tensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<2>::SymTensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<2>::ThirdRankTensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(Dim<3>::Vector&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<3>::Tensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<3>::SymTensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(Dim<3>::ThirdRankTensor&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

    .def("write", (void (Obj::*)(const std::vector<int>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const std::vector<double>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<std::string>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    .def("write", (void (Obj::*)(const std::vector<Dim<1>::Vector>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<1>::Tensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<1>::SymTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<1>::ThirdRankTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    .def("write", (void (Obj::*)(const std::vector<Dim<2>::Vector>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<2>::Tensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<2>::SymTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<2>::ThirdRankTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    .def("write", (void (Obj::*)(const std::vector<Dim<3>::Vector>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<3>::Tensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<3>::SymTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)
    .def("write", (void (Obj::*)(const std::vector<Dim<3>::ThirdRankTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(std::vector<int>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(std::vector<double>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<std::string>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(std::vector<Dim<1>::Vector>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<1>::Tensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<1>::SymTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<1>::ThirdRankTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(std::vector<Dim<2>::Vector>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<2>::Tensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<2>::SymTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<2>::ThirdRankTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

    .def("read", (void (Obj::*)(std::vector<Dim<3>::Vector>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<3>::Tensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<3>::SymTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)
    .def("read", (void (Obj::*)(std::vector<Dim<3>::ThirdRankTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a)

#ifdef SPHERAL1D
    .def("write", (void (Obj::*)(const Field<Dim<1>, Dim<1>::Scalar>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<1>, Dim<1>::Vector>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<1>, Dim<1>::Tensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<1>, Dim<1>::SymTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<1>, Dim<1>::ThirdRankTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<1>, int>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 

    .def("read", (void (Obj::*)(Field<Dim<1>, Dim<1>::Scalar>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<1>, Dim<1>::Vector>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<1>, Dim<1>::Tensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<1>, Dim<1>::SymTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<1>, Dim<1>::ThirdRankTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<1>, int>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
#endif
    
#ifdef SPHERAL2D
    .def("write", (void (Obj::*)(const Field<Dim<2>, Dim<2>::Scalar>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<2>, Dim<2>::Vector>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<2>, Dim<2>::Tensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<2>, Dim<2>::SymTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<2>, Dim<2>::ThirdRankTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<2>, int>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 

    .def("read", (void (Obj::*)(Field<Dim<2>, Dim<2>::Scalar>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<2>, Dim<2>::Vector>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<2>, Dim<2>::Tensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<2>, Dim<2>::SymTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<2>, Dim<2>::ThirdRankTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<2>, int>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
#endif
    
#ifdef SPHERAL1D
    .def("write", (void (Obj::*)(const Field<Dim<3>, Dim<3>::Scalar>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<3>, Dim<3>::Vector>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<3>, Dim<3>::Tensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<3>, Dim<3>::SymTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<3>, Dim<3>::ThirdRankTensor>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 
    .def("write", (void (Obj::*)(const Field<Dim<3>, int>&, const std::string)) &Obj::write, "value"_a, "pathName"_a) 

    .def("read", (void (Obj::*)(Field<Dim<3>, Dim<3>::Scalar>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<3>, Dim<3>::Vector>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<3>, Dim<3>::Tensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<3>, Dim<3>::SymTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<3>, Dim<3>::ThirdRankTensor>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
    .def("read", (void (Obj::*)(Field<Dim<3>, int>&, const std::string) const) &Obj::read, "value"_a, "pathName"_a) 
#endif
    
    .def("write_unsigned_int", (void (Obj::*)(const unsigned, const std::string)) &Obj::write_unsigned_int, "value"_a, "pathName"_a)
    .def("write_int", (void (Obj::*)(const int, const std::string)) &Obj::write_int, "value"_a, "pathName"_a)
    .def("write_bool", (void (Obj::*)(const bool, const std::string)) &Obj::write_bool, "value"_a, "pathName"_a)
    .def("write_double", (void (Obj::*)(const double, const std::string)) &Obj::write_double, "value"_a, "pathName"_a)
    .def("write_string", (void (Obj::*)(const std::string, const std::string)) &Obj::write_string, "value"_a, "pathName"_a)

    .def("read_unsigned_int", (unsigned (Obj::*)(const std::string) const) &Obj::read_unsigned_int, "pathName"_a)
    .def("read_int", (int (Obj::*)(const std::string) const) &Obj::read_int, "pathName"_a)
    .def("read_bool", (bool (Obj::*)(const std::string) const) &Obj::read_bool, "pathName"_a)
    .def("read_double", (double (Obj::*)(const std::string) const) &Obj::read_double, "pathName"_a)
    .def("read_string", (std::string (Obj::*)(const std::string) const) &Obj::read_string, "pathName"_a)

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
  py::class_<FileIO, PyAbstractFileIO<FileIO>> fileioPB11(m, "FileIO");
  virtualFileIOBindings<FileIO>(m, fileioPB11);
  fileioPB11

    // Constructors
    .def(py::init<>())
    .def(py::init<std::string, AccessType>(), "filename"_a, "access"_a)

    // Methods
    .def("writeObject", &FileIO::writeObject)
    .def("readObject", &FileIO::readObject)
    ;

  //............................................................................
  // FlatFileIO
  py::class_<FlatFileIO, FileIO, PyConcreteFileIO<FlatFileIO>> flatfileioPB11(m, "FlatFileIO");
  virtualFileIOBindings<FlatFileIO>(m, flatfileioPB11);
  flatfileioPB11

    // Constructors
    .def(py::init<>())
    .def(py::init<std::string, AccessType, FlatFileFormat>(), "filename"_a, "access"_a, "format"_a=FlatFileFormat::ascii)

    // Methods
    .def("findPathName", &FlatFileIO::findPathName)
    .def("beginingOfFile", &FlatFileIO::beginningOfFile)

    // Attributes
    .def_property("precision", &FlatFileIO::precision, &FlatFileIO::setPrecision)
    .def_property_readonly("readyToWrite", &FlatFileIO::readyToWrite)
    .def_property_readonly("readyToRead", &FlatFileIO::readyToRead)
    ;

  //............................................................................
  // SiloFileIO
  py::class_<SiloFileIO, FileIO, PyConcreteFileIO<SiloFileIO>> silofileioPB11(m, "SiloFileIO");
  virtualFileIOBindings<SiloFileIO>(m, silofileioPB11);
  silofileioPB11

    // Constructors
    .def(py::init<>())
    .def(py::init<std::string, AccessType>(), "filename"_a, "access"_a)
    ;

  //............................................................................
  // PyFileIO
  py::class_<PyFileIO, FileIO, PyPyConcreteFileIO<PyFileIO>> pyfileioPB11(m, "PyFileIO");
  virtualFileIOBindings<PyFileIO>(m, pyfileioPB11);
  pyfileioPB11

    // Constructors
    .def(py::init<>())
    .def(py::init<std::string, AccessType>(), "filename"_a, "access"_a)
    ;

  //............................................................................
  // vector2string
  m.def("vector2string", &vector2string<int>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<unsigned>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<uint64_t>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<double>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<std::string>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<1>::Vector>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<2>::Vector>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<3>::Vector>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<1>::Tensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<2>::Tensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<3>::Tensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<1>::SymTensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<2>::SymTensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<3>::SymTensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<1>::ThirdRankTensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<2>::ThirdRankTensor>, "val"_a, "precision"_a=30);
  m.def("vector2string", &vector2string<Dim<3>::ThirdRankTensor>, "val"_a, "precision"_a=30);

  //............................................................................
  // string2vector
  m.def("string2vector_of_int", &string2vector<int>, "val"_a);
  m.def("string2vector_of_unsigned", &string2vector<unsigned>, "val"_a);
  m.def("string2vector_of_ULL", &string2vector<uint64_t>, "val"_a);
  m.def("string2vector_of_double", &string2vector<double>, "val"_a);
  m.def("string2vector_of_string", &string2vector<std::string>, "val"_a);
  m.def("string2vector_of_Vector1d", &string2vector<Dim<1>::Vector>, "val"_a);
  m.def("string2vector_of_Vector2d", &string2vector<Dim<2>::Vector>, "val"_a);
  m.def("string2vector_of_Vector3d", &string2vector<Dim<3>::Vector>, "val"_a);
  m.def("string2vector_of_Tensor1d", &string2vector<Dim<1>::Tensor>, "val"_a);
  m.def("string2vector_of_Tensor2d", &string2vector<Dim<2>::Tensor>, "val"_a);
  m.def("string2vector_of_Tensor3d", &string2vector<Dim<3>::Tensor>, "val"_a);
  m.def("string2vector_of_SymTensor1d", &string2vector<Dim<1>::SymTensor>, "val"_a);
  m.def("string2vector_of_SymTensor2d", &string2vector<Dim<2>::SymTensor>, "val"_a);
  m.def("string2vector_of_SymTensor3d", &string2vector<Dim<3>::SymTensor>, "val"_a);
  m.def("string2vector_of_ThirdRankTensor1d", &string2vector<Dim<1>::ThirdRankTensor>, "val"_a);
  m.def("string2vector_of_ThirdRankTensor2d", &string2vector<Dim<2>::ThirdRankTensor>, "val"_a);
  m.def("string2vector_of_ThirdRankTensor3d", &string2vector<Dim<3>::ThirdRankTensor>, "val"_a);

  return m.ptr();
}
