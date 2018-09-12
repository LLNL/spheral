// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/EquationOfState.hh"
#include "Material/GammaLawGas.hh"
#include "Material/PolytropicEquationOfState.hh"
#include "Material/IsothermalEquationOfState.hh"
#include "Field/Field.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::Material;
using namespace Spheral::FieldSpace;

namespace Spheral {
namespace Material {

//------------------------------------------------------------------------------
// PyAbstractEOS
//------------------------------------------------------------------------------
template<typename Dimension, class EOSBase>
class PyAbstractEOS: public EOSBase {
public:
  using EOSBase::EOSBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  virtual void setPressure(Field<Dimension, Scalar>& pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setPressure,  // name of method
                           pressure, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setTemperature,  // name of method
                           temperature, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setSpecificThermalEnergy,  // name of method
                           specificThermalEnergy, massDensity, temperature        // arguments
                           );
  }

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setSpecificHeat,  // name of method
                           specificHeat, massDensity, temperature        // arguments
                           );
  }

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setSoundSpeed,  // name of method
                           soundSpeed, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setGammaField(Field<Dimension, Scalar>& gammaField,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setGammaField,  // name of method
                           gammaField, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setBulkModulus,  // name of method
                           bulkModulus, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           EOSBase,      // Parent class
                           setEntropy,  // name of method
                           entropy, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual bool valid() const override {
    PYBIND11_OVERLOAD_PURE(bool,         // Return type
                           EOSBase,      // Parent class
                           valid,  // name of method
                           );
  }

};

//------------------------------------------------------------------------------
// PyEOS
//------------------------------------------------------------------------------
template<typename Dimension, class EOSBase>
class PyEOS: public EOSBase {
public:
  using EOSBase::EOSBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  virtual void setPressure(Field<Dimension, Scalar>& pressure,
                           const Field<Dimension, Scalar>& massDensity,
                           const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setPressure,  // name of method
                           pressure, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setTemperature(Field<Dimension, Scalar>& temperature,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setTemperature,  // name of method
                           temperature, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                                        const Field<Dimension, Scalar>& massDensity,
                                        const Field<Dimension, Scalar>& temperature) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setSpecificThermalEnergy,  // name of method
                           specificThermalEnergy, massDensity, temperature        // arguments
                           );
  }

  virtual void setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                               const Field<Dimension, Scalar>& massDensity,
                               const Field<Dimension, Scalar>& temperature) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setSpecificHeat,  // name of method
                           specificHeat, massDensity, temperature        // arguments
                           );
  }

  virtual void setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setSoundSpeed,  // name of method
                           soundSpeed, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setGammaField(Field<Dimension, Scalar>& gammaField,
                             const Field<Dimension, Scalar>& massDensity,
                             const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setGammaField,  // name of method
                           gammaField, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
                              const Field<Dimension, Scalar>& massDensity,
                              const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setBulkModulus,  // name of method
                           bulkModulus, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual void setEntropy(Field<Dimension, Scalar>& entropy,
                          const Field<Dimension, Scalar>& massDensity,
                          const Field<Dimension, Scalar>& specificThermalEnergy) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           EOSBase,      // Parent class
                           setEntropy,  // name of method
                           entropy, massDensity, specificThermalEnergy        // arguments
                           );
  }

  virtual bool valid() const override {
    PYBIND11_OVERLOAD(bool,         // Return type
                           EOSBase,      // Parent class
                           valid,  // name of method
                           );
  }

};

}
}

namespace {  // anonymous

//------------------------------------------------------------------------------
// Common virtual methods of EOS objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void virtualEOSBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  typedef typename Dimension::Scalar Scalar;

  obj

    // Methods
    .def("setPressure", &Obj::setPressure, "pressure"_a, "massDensity"_a, "specificThermalEnergy"_a)
    .def("setTemperature", &Obj::setTemperature, "temperature"_a, "massDensity"_a, "specificThermalEnergy"_a)
    .def("setSpecificThermalEnergy", &Obj::setSpecificThermalEnergy, "specificThermalEnergy"_a, "massDensity"_a, "temperature"_a)
    .def("setSpecificHeat", &Obj::setSpecificHeat, "specificHeat"_a, "massDensity"_a, "temperature"_a)
    .def("setSoundSpeed", &Obj::setSoundSpeed, "soundSpeed"_a, "massDensity"_a, "specificThermalEnergy"_a)
    .def("setGammaField", &Obj::setGammaField, "gammaField"_a, "massDensity"_a, "specificThermalEnergy"_a)
    .def("setBulkModulus", &Obj::setBulkModulus, "bulkModulus"_a, "massDensity"_a, "specificThermalEnergy"_a)
    .def("setEntropy", &Obj::setEntropy, "entropy"_a, "massDensity"_a, "specificThermalEnergy"_a)
    .def("valid", &Obj::valid)
    ;
}

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  using Spheral::NodeSpace::NodeList;
  using Spheral::FieldSpace::Field;

  //............................................................................
  // EquationOfState
  typedef EquationOfState<Dimension> EOS;
  py::class_<EOS, PyAbstractEOS<Dimension, EOS>> eosPB11(m, ("EquationOfState" + suffix).c_str());
  virtualEOSBindings<Dimension, EOS>(m, suffix, eosPB11);
  eosPB11
    
    // Constructors
    .def(py::init<const PhysicalConstants&, double, double, MaterialPressureMinType>(),
         "constants"_a,
         "minimumPressure"_a = -std::numeric_limits<double>::max(),
         "maximumPressure"_a =  std::numeric_limits<double>::max(),
         "minPressureType"_a = MaterialPressureMinType::PressureFloor)

    // Methods
    .def("applyPressureLimits", &EOS::applyPressureLimits)

    // Attributes
    .def_property_readonly("constants", &EOS::constants)
    .def_property("minimumPressure",
                  (double (EOS::*)() const) &EOS::minimumPressure,
                  (void (EOS::*)(const double)) &EOS::minimumPressure)
    .def_property("maximumPressure",
                  (double (EOS::*)() const) &EOS::maximumPressure,
                  (void (EOS::*)(const double)) &EOS::maximumPressure)
    .def_property("minimumPressureType",
                  (MaterialPressureMinType (EOS::*)() const) &EOS::minimumPressureType,
                  (void (EOS::*)(const MaterialPressureMinType)) &EOS::minimumPressureType)
    ;

  //............................................................................
  // GammaLawGas
  typedef GammaLawGas<Dimension> GLG;
  py::class_<GLG, EOS, PyEOS<Dimension, GLG>> gammalawPB11(m, ("GammaLawGas" + suffix).c_str());
  virtualEOSBindings<Dimension, GLG>(m, suffix, gammalawPB11);
  gammalawPB11

    // Constructors
    .def(py::init<double, double, const PhysicalConstants&, double, double, MaterialPressureMinType>(),
         "gamma"_a,
         "mu"_a,
         "constants"_a,
         "minimumPressure"_a = -std::numeric_limits<double>::max(),
         "maximumPressure"_a =  std::numeric_limits<double>::max(),
         "minPressureType"_a = MaterialPressureMinType::PressureFloor)

    // Methods
    .def("pressure", &GLG::pressure, "massDensity"_a, "specificThermalEnergy"_a)
    .def("temperature", &GLG::temperature, "massDensity"_a, "specificThermalEnergy"_a)
    .def("specificThermalEnergy", &GLG::specificThermalEnergy, "massDensity"_a, "temperature"_a)
    .def("specificHeat", &GLG::specificHeat, "massDensity"_a, "temperature"_a)
    .def("soundSpeed", &GLG::soundSpeed, "massDensity"_a, "specificThermalEnergy"_a)
    .def("gamma", &GLG::gamma, "massDensity"_a, "specificThermalEnergy"_a)
    .def("bulkModulus", &GLG::bulkModulus, "massDensity"_a, "specificThermalEnergy"_a)
    .def("entropy", &GLG::entropy, "massDensity"_a, "specificThermalEnergy"_a)

    // Attributes
    .def_property("gamma", &GLG::getGamma, &GLG::setGamma)
    .def_property("mu", &GLG::getMolecularWeight, &GLG::setMolecularWeight)
    ;

  //............................................................................
  // PolytropicEquationOfState
  typedef PolytropicEquationOfState<Dimension> PolyEOS;
  py::class_<PolyEOS, EOS, PyEOS<Dimension, PolyEOS>> polyeosPB11(m, ("PolytropicEquationOfState" + suffix).c_str());
  virtualEOSBindings<Dimension, PolyEOS>(m, suffix, polyeosPB11);
  polyeosPB11

    // Constructors
    .def(py::init<double, double, double, const PhysicalConstants&, double, double, MaterialPressureMinType>(),
         "K"_a,
         "index"_a,
         "mu"_a,
         "constants"_a,
         "minimumPressure"_a = -std::numeric_limits<double>::max(),
         "maximumPressure"_a =  std::numeric_limits<double>::max(),
         "minPressureType"_a = MaterialPressureMinType::PressureFloor)

    // Methods
    .def("pressure", &PolyEOS::pressure, "massDensity"_a, "specificThermalEnergy"_a)
    .def("temperature", &PolyEOS::temperature, "massDensity"_a, "specificThermalEnergy"_a)
    .def("specificThermalEnergy", &PolyEOS::specificThermalEnergy, "massDensity"_a, "temperature"_a)
    .def("specificHeat", &PolyEOS::specificHeat, "massDensity"_a, "temperature"_a)
    .def("soundSpeed", &PolyEOS::soundSpeed, "massDensity"_a, "specificThermalEnergy"_a)
    .def("gamma", (Scalar (PolyEOS::*)(const Scalar, const Scalar) const) &PolyEOS::gamma, "massDensity"_a, "specificThermalEnergy"_a)
    .def("bulkModulus", &PolyEOS::bulkModulus, "massDensity"_a, "specificThermalEnergy"_a)
    .def("entropy", &PolyEOS::entropy, "massDensity"_a, "specificThermalEnergy"_a)

    // Attributes
    .def_property_readonly("polytropicConstant", &PolyEOS::polytropicConstant)
    .def_property_readonly("polytropicIndex", &PolyEOS::polytropicIndex)
    .def_property_readonly("gamma_", (double (PolyEOS::*)() const) &PolyEOS::gamma)
    .def_property_readonly("molecularWeight", &PolyEOS::molecularWeight)
    .def_property("externalPressure", &PolyEOS::externalPressure, &PolyEOS::setExternalPressure)
    ;

  //............................................................................
  // IsothermalEquationOfState
  typedef IsothermalEquationOfState<Dimension> IsoEOS;
  py::class_<IsoEOS, EOS, PyEOS<Dimension, IsoEOS>> isoeosPB11(m, ("IsothermalEquationOfState" + suffix).c_str());
  virtualEOSBindings<Dimension, IsoEOS>(m, suffix, isoeosPB11);
  isoeosPB11

    // Constructors
    .def(py::init<double, double, const PhysicalConstants&, double, double, MaterialPressureMinType>(),
         "K"_a,
         "mu"_a,
         "constants"_a,
         "minimumPressure"_a = -std::numeric_limits<double>::max(),
         "maximumPressure"_a =  std::numeric_limits<double>::max(),
         "minPressureType"_a = MaterialPressureMinType::PressureFloor)

    // Methods
    .def("pressure", &IsoEOS::pressure, "massDensity"_a, "specificThermalEnergy"_a)
    .def("temperature", &IsoEOS::temperature, "massDensity"_a, "specificThermalEnergy"_a)
    .def("specificThermalEnergy", &IsoEOS::specificThermalEnergy, "massDensity"_a, "temperature"_a)
    .def("specificHeat", &IsoEOS::specificHeat, "massDensity"_a, "temperature"_a)
    .def("soundSpeed", &IsoEOS::soundSpeed, "massDensity"_a, "specificThermalEnergy"_a)
    .def("gamma", (Scalar (IsoEOS::*)(const Scalar, const Scalar) const) &IsoEOS::gamma, "massDensity"_a, "specificThermalEnergy"_a)
    .def("bulkModulus", &IsoEOS::bulkModulus, "massDensity"_a, "specificThermalEnergy"_a)
    .def("entropy", &IsoEOS::entropy, "massDensity"_a, "specificThermalEnergy"_a)

    // Attributes
    .def_property_readonly("K", &IsoEOS::K)
    .def_property_readonly("molecularWeight", &IsoEOS::molecularWeight)
    .def_property("externalPressure", &IsoEOS::externalPressure, &IsoEOS::setExternalPressure)
    ;

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralMaterial) {
  py::module m("SpheralMaterial", "Spheral Material module.");

  //............................................................................
  // MaterialPressureMinType
  py::enum_<Spheral::MaterialPressureMinType>(m, "MaterialPressureMinType")
    .value("PressureFloor", Spheral::MaterialPressureMinType::PressureFloor)
    .value("ZeroPressure",  Spheral::MaterialPressureMinType::ZeroPressure)
    .export_values();

  //............................................................................
  // PhysicalConstants
  py::class_<PhysicalConstants>(m, "PhysicalConstants")

    // Constructors
    .def(py::init<double, double, double>(), "unitLm"_a, "unitMkg"_a, "unitTsec"_a)

    // Attributes
    .def_property_readonly("unitLengthMeters", &PhysicalConstants::unitLengthMeters)
    .def_property_readonly("unitMassKg", &PhysicalConstants::unitMassKg)
    .def_property_readonly("unitTimeSec", &PhysicalConstants::unitTimeSec)
    .def_property_readonly("unitMassDensity", &PhysicalConstants::unitMassDensity)
    .def_property_readonly("protonMass", &PhysicalConstants::protonMass)
    .def_property_readonly("electronMass", &PhysicalConstants::electronMass)
    .def_property_readonly("electronCharge", &PhysicalConstants::electronCharge)
    .def_property_readonly("G", &PhysicalConstants::G)
    .def_property_readonly("c", &PhysicalConstants::c)
    .def_property_readonly("kB", &PhysicalConstants::kB)
    .def_property_readonly("Navogadro", &PhysicalConstants::Navogadro)
    .def_property_readonly("molarGasConstant", &PhysicalConstants::molarGasConstant)
    .def_property_readonly("kelvinsToEnergyPerMole", &PhysicalConstants::kelvinsToEnergyPerMole)
    .def_property_readonly("stefanBoltzmannConstant", &PhysicalConstants::stefanBoltzmannConstant)
    ;

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif
#ifdef SPHERAL2D
  dimensionBindings<Spheral::Dim<2>>(m, "2d");
#endif
#ifdef SPHERAL3D
  dimensionBindings<Spheral::Dim<3>>(m, "3d");
#endif

  return m.ptr();
}
