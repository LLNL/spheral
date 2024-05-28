//---------------------------------Spheral++----------------------------------//
// SubPointPressureHourglassControl
//
// Impose additional forces on each point using subdivisions of the Voronoi
// control volume to constrain the unphysical degrees of freedom in our hydro
// discretization and avoid spurious so-called houglass modes.
//----------------------------------------------------------------------------//
#include "VoronoiCells/SubPointPressureHourglassControl.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/CellFaceFlag.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Utilities/Timer.hh"

#include <limits>
#include <algorithm>

namespace Spheral {

using std::vector;

namespace {  // anonymous

template<typename Vector>
inline
double
limiter(const Vector& xi,
        const Vector& xn,
        const double  rhoi,
        const Vector& gradRhoi) {
  const auto delta = gradRhoi.dot(xn - xi);
  return (rhoi + delta < 0.0 ?
          abs(rhoi/delta) :
          1.0);
}

//------------------------------------------------------------------------------
// Center of mass (1D)
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
centerOfMass(const Dim<1>::FacetedVolume& celli,
             const Dim<1>::Vector& xi,
             const Dim<1>::Scalar rhoi,
             Dim<1>::Vector gradRhoi) {
  using Vector = Dim<1>::Vector;
  const auto& x1 = celli.xmin();
  const auto& x2 = celli.xmax();
  gradRhoi *= limiter(xi, x1, rhoi, gradRhoi);
  gradRhoi *= limiter(xi, x2, rhoi, gradRhoi);
  CHECK(rhoi + gradRhoi.dot(x1 - xi) >= 0.0);
  CHECK(rhoi + gradRhoi.dot(x2 - xi) >= 0.0);
  const auto c = xi.x();
  const auto ab = x1.x() + x2.x() - 2.0*c;
  const auto m = gradRhoi.x();
  const Vector result(c + ab*(2.0*m*ab + 3.0*rhoi)/(3.0*(m*ab + 2.0*rhoi)));
  ENSURE2(result.x() >= x1.x() and result.x() <= x2.x(), result << " not in [" << x1.x() << " " << x2.x() << "] : " << gradRhoi << " " << (rhoi + gradRhoi.dot(x1 - xi)) << " " << (rhoi + gradRhoi.dot(x2 - xi)));
  return result;
}

//------------------------------------------------------------------------------
// Center of mass (2D)
//------------------------------------------------------------------------------
inline
Dim<2>::Vector
centerOfMass(const Dim<2>::FacetedVolume& celli,
             const Dim<2>::Vector& xi,
             const Dim<2>::Scalar rhoi,
             Dim<2>::Vector gradRhoi) {
  using Vector = Dim<2>::Vector;
  return Vector();
}

//------------------------------------------------------------------------------
// Center of mass (3D)
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
centerOfMass(const Dim<3>::FacetedVolume& celli,
             const Dim<3>::Vector& xi,
             const Dim<3>::Scalar rhoi,
             Dim<3>::Vector gradRhoi) {
  using Vector = Dim<3>::Vector;
  return Vector();
}

//------------------------------------------------------------------------------
// Compute the internal acceleration (1D)
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
subCellAcceleration(const Dim<1>::FacetedVolume& celli,
                    const Dim<1>::Vector& xi,
                    const Dim<1>::Scalar  Pi,
                    const Dim<1>::Scalar  rhoi,
                    const Dim<1>::Vector& gradRhoi) {
  using Vector = Dim<1>::Vector;
  const auto comi = centerOfMass(celli, xi, rhoi, gradRhoi); // celli.centroid();

  // Define a function to increment the acceleration for each subcell
  auto asub = [&](const Vector& vert) -> Vector {
                const auto dA = (comi - vert).unitVector();
                const auto Psub = abs(Pi * (vert.x() - comi.x())/(vert.x() - xi.x()));
                return Psub*dA;
              };

  // Now we can sum up finite volume contribution to the acceleration for each subvolume
  const auto Vi = celli.volume();
  CHECK(Vi > 0.0);
  return (asub(celli.xmin()) + asub(celli.xmax()))/(rhoi*Vi);
}

//------------------------------------------------------------------------------
// Compute the internal acceleration (2D)
//------------------------------------------------------------------------------
inline
Dim<2>::Vector
subCellAcceleration(const Dim<2>::FacetedVolume& celli,
                    const Dim<2>::Vector& xi,
                    const Dim<2>::Scalar  Pi,
                    const Dim<2>::Scalar  rhoi,
                    const Dim<2>::Vector& gradRhoi) {
  using Vector = Dim<2>::Vector;
  const auto comi = celli.centroid();

  // Define a function to increment the acceleration for each subcell
  auto asub = [&](const Vector& v1, const Vector& v2) -> Vector {
                const auto v12 = v2 - v1;
                const Vector dA(-v12.y(), v12.x());
                const auto Psub = abs(Pi * ((v1 - comi).cross(v2 - comi)).z()*safeInv(((v1 - xi).cross(v2 - xi)).z()));
                return Psub*dA;
              };

  // Now we can sum up finite volume contribution to the acceleration for each subvolume.
  Vector result;
  const auto& facets = celli.facets();
  for (auto& f: facets) result += asub(f.point1(), f.point2());
  const auto Vi = celli.volume();
  CHECK(Vi > 0.0);
  result /= rhoi*Vi;
  return result;
}

//------------------------------------------------------------------------------
// Compute the internal acceleration (3D)
//------------------------------------------------------------------------------
inline
Dim<3>::Vector
subCellAcceleration(const Dim<3>::FacetedVolume& celli,
                    const Dim<3>::Vector& xi,
                    const Dim<3>::Scalar  Pi,
                    const Dim<3>::Scalar  rhoi,
                    const Dim<3>::Vector& gradRhoi) {
  using Vector = Dim<3>::Vector;
  return Vector();
}

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
SubPointPressureHourglassControl<Dimension>::
SubPointPressureHourglassControl(const Scalar fHG):
  mfHG(fHG) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SubPointPressureHourglassControl<Dimension>::
~SubPointPressureHourglassControl() {
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
}

//------------------------------------------------------------------------------
// No derivatives to register
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// No time step vote
//------------------------------------------------------------------------------
template<typename Dimension>
typename SubPointPressureHourglassControl<Dimension>::TimeStepType
SubPointPressureHourglassControl<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {
  return std::make_pair(std::numeric_limits<double>::max(), std::string("SubPointPressureHourglassControl: no vote"));
}

//------------------------------------------------------------------------------
// Add our terms to the hydro derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
SubPointPressureHourglassControl<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_BEGIN("SubPointHGevalDerivs");

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  numNodeLists = nodeLists.size();
  const auto  npairs = pairs.size();
  const auto  nint = dataBase.numInternalNodes();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto vel = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto P = state.fields(HydroFieldNames::pressure, 0.0);
  const auto cells = state.template fields<FacetedVolume>(HydroFieldNames::cells);
  const auto cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, CellFaceFlag());
  const auto gradRho = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  CHECK(mass.size() == numNodeLists);
  CHECK(pos.size() == numNodeLists);
  CHECK(rho.size() == numNodeLists);
  CHECK(P.size() == numNodeLists);
  CHECK(cells.size() == numNodeLists);
  CHECK(cellFaceFlags.size() == numNodeLists);
  CHECK(gradRho.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto& pairAccelerations = derivs.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  const auto compatibleEnergy = pairAccelerations.size() > 0u;
  if (compatibleEnergy) {
    CHECK(pairAccelerations.size() == npairs or pairAccelerations.size() == npairs + nint);
    if (pairAccelerations.size() == npairs) {
      pairAccelerations.resize(npairs + nint);
      std::fill(pairAccelerations.begin() + npairs, pairAccelerations.end(), Vector::zero);
    }
  }

  // Walk the points
  auto offset = npairs;
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto n = mass[k]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto& xi = pos(k,i);
      const auto& vi = vel(k,i);
      const auto& celli = cells(k,i);
      // const auto  mi = mass(k,i);
      const auto  Pi = P(k,i);
      const auto  rhoi = rho(k,i);
      const auto& gradRhoi = gradRho(k,i);
      const auto  deltaDvDti = mfHG * subCellAcceleration(celli, xi, Pi, rhoi, gradRhoi);
      DvDt(k,i) += deltaDvDti;
      DepsDt(k,i) -= vi.dot(deltaDvDti);
      if (compatibleEnergy) pairAccelerations[offset + i] += deltaDvDti;
    }
    offset += n;
  }

  TIME_END("SubPointHGevalDerivs");
}

} // end namespace Spheral
