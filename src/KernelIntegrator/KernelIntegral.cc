//---------------------------------Spheral++----------------------------------//
// KernelIntegral
//
// Represents an integral to be integrated by Integrator
//----------------------------------------------------------------------------//
#include "KernelIntegral.hh"

#include <cmath>
#include <limits>

#include "Hydro/HydroFieldNames.hh"

namespace Spheral {

constexpr auto epsilon = std::numeric_limits<double>::epsilon() * 10;

//------------------------------------------------------------------------------
// LinearKernel
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearKernel<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    this->mValues[locali] += kid.weight * coeff * kid.values[i];
  }
}

  
//------------------------------------------------------------------------------
// LinearGrad
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearGrad<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    this->mValues[locali] += kid.weight * coeff * kid.dvalues[i];
  }
}

//------------------------------------------------------------------------------
// LinearKernelVector
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearKernelVector<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    this->mValues[locali] += kid.weight * coeff * kid.values[i];
  }
}

//------------------------------------------------------------------------------
// LinearKernelStdVector
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearKernelStdVector<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  CHECK(coeff.size() == mSize);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    auto& vals = this->mValues[locali];
    CHECK(vals.size() == mSize);
    for (auto k = 0u; k < mSize; ++k) {
      vals[k] += kid.weight * coeff[k] * kid.values[i];
    }
  }
}

//------------------------------------------------------------------------------
// LinearGradStdVector
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearGradStdVector<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  CHECK(coeff.size() == mSize);
  const auto numIndices = kid.indices.size();
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    auto& vals = this->mValues[locali];
    CHECK(vals.size() == mSize);
    for (auto k = 0u; k < mSize; ++k) {
      vals[k] += kid.weight * coeff[k] * kid.dvalues[i];
    }
  }
}

//------------------------------------------------------------------------------
// BilinearKernelKernel
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearKernelKernel<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (std::abs(kid.values[i]) > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * kid.values[i] * kid.values[j];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearGradKernel
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearGradKernel<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (kid.dvalues[i].magnitude2() > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * kid.dvalues[i] * kid.values[j];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearKernelGrad
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearKernelGrad<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (std::abs(kid.values[i]) > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * kid.values[i] * kid.dvalues[j];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearGradProdGrad
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearGradProdGrad<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (kid.dvalues[i].magnitude2() > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * kid.dvalues[i].dyad(kid.dvalues[j]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearGradDotGrad
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearGradDotGrad<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (kid.dvalues[i].magnitude2() > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * kid.dvalues[i].dot(kid.dvalues[j]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearSurfaceNormalKernelKernelFromGrad
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearSurfaceNormalKernelKernelFromGrad<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (kid.dvalues[i].magnitude2() + std::abs(kid.values[i]) > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * (kid.dvalues[i] * kid.values[j] + kid.values[i] * kid.dvalues[j]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// LinearSurfaceKernel
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSurfaceKernel<Dimension>::
addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.surfaceIndex.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    this->mValues[locali] += kid.weight * coeff * kid.values[i];
  }
}

//------------------------------------------------------------------------------
// LinearSurfaceNormalKernel
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSurfaceNormalKernel<Dimension>::
addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.surfaceIndex.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    auto& vals = this->mValues[locali];
    const auto s = kid.surfaceIndex[i];
    if (s != FlatConnectivity<Dimension>::NoConnectivity) {
      CHECK(size_t(s) < vals.size());
      vals[s] += kid.weight * coeff * kid.normal * kid.values[i];
    }
  }
}

//------------------------------------------------------------------------------
// LinearSurfaceNormalKernelStdVector
//------------------------------------------------------------------------------
template<typename Dimension>
void
LinearSurfaceNormalKernelStdVector<Dimension>::
addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  CHECK(coeff.size() == mSize);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.surfaceIndex.size() == numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    auto& vals = this->mValues[locali];
    const auto s = kid.surfaceIndex[i];
    if (s != FlatConnectivity<Dimension>::NoConnectivity) {
      CHECK(size_t(s) < vals.size());
      CHECK(vals[s].size() == mSize);
      for (auto k = 0u; k < mSize; ++k) {
        vals[s][k] += kid.weight * coeff[k] * kid.normal * kid.values[i];
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearSurfaceKernelKernel
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearSurfaceKernelKernel<Dimension>::
addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (std::abs(kid.values[i]) > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * kid.values[i] * kid.values[j];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearSurfaceNormalKernelKernel
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearSurfaceNormalKernelKernel<Dimension>::
addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.surfaceIndex.size() == numIndices);
  CHECK(kid.numSurfaces.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (std::abs(kid.values[i]) > epsilon) {
      auto& vals = this->mValues[locali];
      const auto s = kid.surfaceIndex[i];
      if (s != FlatConnectivity<Dimension>::NoConnectivity) {
        const auto numSurfaces = kid.numSurfaces[i];
        CHECK(s < numSurfaces);
        for (auto j = 0u; j < numIndices; ++j) {
          const auto k = kid.localIndex[j + numIndices * i];
          if (k != FlatConnectivity<Dimension>::NoConnectivity) {
            const auto jslocal = s + numSurfaces * k;
            CHECK(size_t(jslocal) < vals.size());
            vals[jslocal] += kid.weight * coeff * kid.values[i] * kid.values[j] * kid.normal;
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// BilinearSurfaceNormalKernelDotGrad
//------------------------------------------------------------------------------
template<typename Dimension>
void
BilinearSurfaceNormalKernelDotGrad<Dimension>::
addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto numIndices = kid.indices.size();
  CHECK(kid.values.size() == numIndices);
  CHECK(kid.dvalues.size() == numIndices);
  CHECK(kid.indices.size() == numIndices);
  CHECK(kid.localIndex.size() == numIndices * numIndices);
  for (auto i = 0u; i < numIndices; ++i) {
    const size_t locali = kid.indices[i];
    CHECK(locali < this->mValues.size());
    if (kid.dvalues[i].magnitude2() > epsilon) {
      auto& vals = this->mValues[locali];
      for (auto j = 0u; j < numIndices; ++j) {
        const auto k = kid.localIndex[j + numIndices * i];
        if (k != FlatConnectivity<Dimension>::NoConnectivity) {
      
          CHECK(size_t(k) < vals.size());
          vals[k] += kid.weight * coeff * kid.values[i] * kid.dvalues[j].dot(kid.normal);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// CellCoefficient
//------------------------------------------------------------------------------
template<typename Dimension>
void
CellCoefficient<Dimension>::
addToIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const size_t locali = kid.index0;
  CHECK(locali < this->mValues.size());
  this->mValues[locali] += kid.weight * coeff;
}

//------------------------------------------------------------------------------
// SurfaceNormalCoefficient
//------------------------------------------------------------------------------
template<typename Dimension>
void
SurfaceNormalCoefficient<Dimension>::
addToSurfaceIntegral(const KernelIntegrationData<Dimension>& kid) {
  const auto coeff = this->mCoefficient->evaluateCoefficient(kid);
  const auto s = kid.surfaceIndex0;
  CHECK(s != FlatConnectivity<Dimension>::NoConnectivity);
  const size_t locali = kid.index0;
  CHECK(locali < this->mValues.size());
  auto& vals = this->mValues[locali];
  CHECK(size_t(s) < vals.size());
  vals[s] += kid.weight * coeff * kid.normal;
}

} // end namespace Spheral
