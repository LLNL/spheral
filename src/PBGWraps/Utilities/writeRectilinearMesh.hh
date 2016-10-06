#ifndef __PBGWRAPS_WRITERECTILINEARMESH__
#define __PBGWRAPS_WRITERECTILINEARMESH__

//---------------------------------Spheral++----------------------------------//
// writeRectilinearMesh
//
// Expose the third party visit dumping routines provided by Hank Childs.
//----------------------------------------------------------------------------//
#include <vector>
#include "Python.h"

#include "visit_writer.h"

#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

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
inline
void
writeRectilinearMesh(const std::string& fileName,
                     const bool binaryFile,
                     const std::vector<int>& dimensions,
                     const std::vector<std::vector<double> >& coords,
                     const std::vector<std::string>& scalarNames,
                     const std::vector<std::string>& vectorNames,
                     const std::vector<std::string>& tensorNames,
                     const std::vector<std::string>& symTensorNames,
                     const std::vector<std::vector<double> >& sampledScalars,
                     const std::vector<std::vector<typename Dimension::Vector> >& sampledVectors,
                     const std::vector<std::vector<typename Dimension::Tensor> >& sampledTensors,
                     const std::vector<std::vector<typename Dimension::SymTensor> >& sampledSymTensors) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Pre-conditions.
  VERIFY(dimensions.size() == 3);
  const int nx = dimensions[0];
  const int ny = dimensions[1];
  const int nz = dimensions[2];
  const int n = nx*ny*nz;
  VERIFY(nx > 0);
  VERIFY(ny > 0);
  VERIFY(nz > 0);
  VERIFY(n > 0);

  VERIFY(coords.size() == 3);
  VERIFY(coords[0].size() == nx + 1);
  VERIFY(coords[1].size() == ny + 1);
  VERIFY(coords[2].size() == nz + 1);

  VERIFY(scalarNames.size() == sampledScalars.size());
  VERIFY(vectorNames.size() == sampledVectors.size());
  VERIFY(tensorNames.size() == sampledTensors.size());
  VERIFY(symTensorNames.size() == sampledSymTensors.size());

  const int nvars = (sampledScalars.size() + 
                     sampledVectors.size() +
                     sampledTensors.size() +
                     sampledSymTensors.size());

  for (int k = 0; k != sampledScalars.size(); ++k) VERIFY2(sampledScalars[k].size() == n, "Wrong size:  " << sampledScalars[k].size() << " " << n);
  for (int k = 0; k != sampledVectors.size(); ++k) VERIFY2(sampledVectors[k].size() == n,  "Wrong size:  " << sampledVectors[k].size() << " " << n);
  for (int k = 0; k != sampledTensors.size(); ++k) VERIFY2(sampledTensors[k].size() == n,  "Wrong size:  " << sampledTensors[k].size() << " " << n);
  for (int k = 0; k != sampledSymTensors.size(); ++k) VERIFY2(sampledSymTensors[k].size() == n, "Wrong size:  " << sampledSymTensors[k].size() << " " << n);

  // Build the visit representation of the dimensions and coordinates.
  std::vector<int> dims(3);
  dims[0] = nx + 1;
  dims[1] = ny + 1;
  dims[2] = nz + 1;

  std::vector<float> x(nx + 1);
  std::vector<float> y(ny + 1);
  std::vector<float> z(nz + 1);
  for (int i = 0; i != nx + 1; ++i) x[i] = float(coords[0][i]);
  for (int i = 0; i != ny + 1; ++i) y[i] = float(coords[1][i]);
  for (int i = 0; i != nz + 1; ++i) z[i] = float(coords[2][i]);

  // Now build the visit representations for the field values.
  std::vector<const char*> varnames;
  std::vector<int> vardim;
  float** vars = new float* [nvars];

  // Extract the scalar fields.
  int varindex = 0;
  {
    for (int k = 0; k != sampledScalars.size(); ++k, ++varindex) {
      varnames.push_back(scalarNames[k].c_str());
      vars[varindex] = new float[n];
      vardim.push_back(1);
      for (int i = 0; i != n; ++i) vars[varindex][i] = float(sampledScalars[k][i]);
    }
  }

  // Extract the vector fields.
  {
    for (int k = 0; k != sampledVectors.size(); ++k, ++varindex) {
      varnames.push_back(vectorNames[k].c_str());
      vars[varindex] = new float[3*n];
      vardim.push_back(3);
      int j = 0;
      for (int i = 0; i != n; ++i) {
        const Dim<3>::Vector vi = convertToVector3d(sampledVectors[k][i]);
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

  // Extract the tensor fields.
  {
    for (int k = 0; k != sampledTensors.size(); ++k, ++varindex) {
      varnames.push_back(tensorNames[k].c_str());
      vars[varindex] = new float[9*n];
      vardim.push_back(9);
      int j = 0;
      for (int i = 0; i != n; ++i) {
        const Dim<3>::Tensor ti = convertToTensor3d(sampledTensors[k][i]);
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

  // Extract the symmetric tensor fields.
  {
    for (int k = 0; k != sampledSymTensors.size(); ++k, ++varindex) {
      varnames.push_back(symTensorNames[k].c_str());
      vars[varindex] = new float[9*n];
      vardim.push_back(9);
      int j = 0;
      for (int i = 0; i != n; ++i) {
        const Dim<3>::Tensor ti = convertToTensor3d(sampledSymTensors[k][i]);
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

#endif
