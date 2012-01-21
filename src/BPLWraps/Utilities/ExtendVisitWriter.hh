//---------------------------------Spheral++----------------------------------//
// ExtendVisitWriter
//
// Expose the third party visit dumping routines provided by Hank Childs.
//----------------------------------------------------------------------------//
#include <vector>

#include "boost/python.hpp"
#include "boost/python/detail/api_placeholder.hpp" // This is where len lives!!!

#include "visit_writer.h"

#include "Geometry/Dimension.hh"
#include "DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Convert the input Vector to a Vector3d.
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
convertToVector3d(const Dim<2>::Vector& x) {
  return Dim<3>::Vector(x.x(), x.y(), 0.0);
}

inline
Dim<3>::Vector
convertToVector3d(const Dim<3>::Vector& x) {
  return x;
}

//------------------------------------------------------------------------------
// Convert the input tensor to a Tensor3d.
//------------------------------------------------------------------------------
inline
Dim<3>::Tensor
convertToTensor3d(const Dim<2>::Tensor& x) {
  return Dim<3>::Tensor(x.xx(), x.xy(), 0.0,
                        x.yx(), x.yy(), 0.0,
                        0.0,    0.0,    0.0);
}

inline
Dim<3>::Tensor
convertToTensor3d(const Dim<2>::SymTensor& x) {
  return Dim<3>::Tensor(x.xx(), x.xy(), 0.0,
                        x.yx(), x.yy(), 0.0,
                        0.0,    0.0,    0.0);
}

inline
Dim<3>::Tensor
convertToTensor3d(const Dim<3>::Tensor& x) {
  return x;
}

inline
Dim<3>::Tensor
convertToTensor3d(const Dim<3>::SymTensor& x) {
  return Dim<3>::Tensor(x.xx(), x.xy(), x.xz(),
                        x.yx(), x.yy(), x.yz(),
                        x.zx(), x.zy(), x.zz());
}

//------------------------------------------------------------------------------
// General templated method.
//------------------------------------------------------------------------------
template<typename Dimension>
void
writeRectilinearMesh(const std::string& fileName,
                     const bool binaryFile,
                     const boost::python::tuple& dimensions,
                     const boost::python::tuple& coords,
                     const boost::python::tuple& sampledFields) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  VERIFY(boost::python::len(dimensions) == 3);
  const int nx = boost::python::extract<int>(dimensions[0]);
  const int ny = boost::python::extract<int>(dimensions[1]);
  const int nz = boost::python::extract<int>(dimensions[2]);
  const int n = nx*ny*nz;
  VERIFY(nx > 0);
  VERIFY(ny > 0);
  VERIFY(nz > 0);
  VERIFY(n > 0);

  VERIFY(boost::python::len(coords) == 3);
  VERIFY(boost::python::len(coords[0]) == nx + 1);
  VERIFY(boost::python::len(coords[1]) == ny + 1);
  VERIFY(boost::python::len(coords[2]) == nz + 1);

  VERIFY(boost::python::len(sampledFields) == 4);
  int nvars = 0;
  for (int i = 0; i != 4; ++i) {
    const int ni = boost::python::len(sampledFields[i]);
    nvars += ni;
    for (int j = 0; j != ni; ++j) 
      VERIFY(boost::python::len(sampledFields[i][j]) == 2);
  }

  // Build the visit representation of the dimensions and coordinates.
  std::vector<int> dims(3);
  dims[0] = nx + 1;
  dims[1] = ny + 1;
  dims[2] = nz + 1;

  std::vector<float> x(nx + 1);
  std::vector<float> y(ny + 1);
  std::vector<float> z(nz + 1);
  for (int i = 0; i != nx + 1; ++i)
    x[i] = boost::python::extract<float>(coords[0][i]);
  for (int i = 0; i != ny + 1; ++i)
    y[i] = boost::python::extract<float>(coords[1][i]);
  for (int i = 0; i != nz + 1; ++i)
    z[i] = boost::python::extract<float>(coords[2][i]);

  // Now build the visit representations for the field values.
  std::vector<const char*> varnames;
  std::vector<int> vardim;
  float** vars = new float* [nvars];

  // Extract the scalar fields.
  int varindex = 0;
  {
    const int ni = boost::python::len(sampledFields[0]);
    for (int k = 0; k != ni; ++k, ++varindex) {
      CHECK(varindex < nvars);
      vars[varindex] = new float[n];
      const char* name = boost::python::extract<char*>(sampledFields[0][k][0]);
      varnames.push_back(name);
      vardim.push_back(1);
      const std::vector<double>& var = boost::python::extract<std::vector<double>&>(sampledFields[0][k][1]);
      CHECK(var.size() == n);
      for (int i = 0; i != n; ++i) vars[varindex][i] = float(var[i]);
    }
  }

  // Extract the vector fields.
  {
    const int ni = boost::python::len(sampledFields[1]);
    for (int k = 0; k != ni; ++k, ++varindex) {
      CHECK(varindex < nvars);
      vars[varindex] = new float[3*n];
      const char* name = boost::python::extract<char*>(sampledFields[1][k][0]);
      varnames.push_back(name);
      vardim.push_back(3);
      const std::vector<Vector>& var = boost::python::extract<std::vector<Vector>&>(sampledFields[1][k][1]);
      CHECK(var.size() == n);
      int j = 0;
      for (int i = 0; i != n; ++i) {
        const Dim<3>::Vector vi = convertToVector3d(var[i]);
        for (Dim<3>::Vector::const_iterator itr = vi.begin();
             itr != vi.end();
             ++itr, ++j) {
          CHECK(j < 3*n);
          vars[varindex][j] = float(*itr);
        }
      }
      CHECK(j == 3*n);
    }
  }

  // Extract the Tensor fields.
  {
    const int ni = boost::python::len(sampledFields[2]);
    for (int k = 0; k != ni; ++k, ++varindex) {
      CHECK(varindex < nvars);
      vars[varindex] = new float[9*n];
      const char* name = boost::python::extract<char*>(sampledFields[2][k][0]);
      varnames.push_back(name);
      vardim.push_back(9);
      const std::vector<Tensor>& var = boost::python::extract<std::vector<Tensor>&>(sampledFields[2][k][1]);
      CHECK(var.size() == n);
      int j = 0;
      for (int i = 0; i != n; ++i) {
        const Dim<3>::Tensor ti = convertToTensor3d(var[i]);
        for (Dim<3>::Tensor::const_iterator itr = ti.begin();
             itr != ti.end();
             ++itr, ++j) {
          CHECK(j < 9*n);
          vars[varindex][j] = float(*itr);
        }
      }
      CHECK(j == 9*n);
    }
  }

  // Extract the SymTensor fields.
  {
    const int ni = boost::python::len(sampledFields[3]);
    for (int k = 0; k != ni; ++k, ++varindex) {
      CHECK(varindex < nvars);
      vars[varindex] = new float[9*n];
      const char* name = boost::python::extract<char*>(sampledFields[3][k][0]);
      varnames.push_back(name);
      vardim.push_back(9);
      const std::vector<SymTensor>& var = boost::python::extract<std::vector<SymTensor>&>(sampledFields[3][k][1]);
      CHECK(var.size() == n);
      int j = 0;
      for (int i = 0; i != n; ++i) {
        const Dim<3>::Tensor ti = convertToTensor3d(var[i]);
        for (Dim<3>::Tensor::const_iterator itr = ti.begin();
             itr != ti.end();
             ++itr, ++j) {
          CHECK(j < 9*n);
          vars[varindex][j] = float(*itr);
        }
      }
      CHECK(j == 9*n);
    }
  }

  CHECK(varindex == nvars);

  // We list all variables as cell centered.
  std::vector<int> centering(nvars, 0);

  // Binary flag.
  const int useBinary = binaryFile ? 1: 0;

  // We now have all the data in visit friendly format, so dump it out.
  CHECK(dims.size() == 3);
  CHECK(x.size() == nx + 1);
  CHECK(y.size() == ny + 1);
  CHECK(z.size() == nz + 1);
  CHECK(vardim.size() == nvars);
  CHECK(centering.size() == nvars);
  write_rectilinear_mesh(fileName.c_str(),
                         useBinary,
                         &(*dims.begin()),
                         &(*x.begin()),
                         &(*y.begin()),
                         &(*z.begin()),
                         nvars,
                         &(*vardim.begin()),
                         &(*centering.begin()),
                         &(*varnames.begin()),
                         vars);

  // Sigh.  Now clean up the space we allocated for vars.
  for (int i = 0; i != nvars; ++i) delete [] vars[i];
  delete [] vars;

}

}
