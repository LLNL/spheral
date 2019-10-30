set(Hydro_inst
  EntropyPolicy
  GammaPolicy
  NonSymmetricSpecificThermalEnergyPolicy
  PositionPolicy
  PressurePolicy
  SecondMomentHourglassControl
  SoundSpeedPolicy
  SpecificFromTotalThermalEnergyPolicy
  SpecificThermalEnergyPolicy
  SumVoronoiMassDensityPolicy
  ThirdMomentHourglassControl
  VolumePolicy
  VoronoiMassDensityPolicy
  )

set(Hydro_sources
  HydroFieldNames.cc
  #TotalHydro.cc
  #VoronoiHourglassControl.cc
  #MashCorrectionPolicy.cc
  )

instantiate(Hydro_inst Hydro_sources)

set(Hydro_headers
  EntropyPolicy.hh
  entropyWeightingFunction.hh
  GammaPolicy.hh
  HydroFieldNames.hh
  HydroInline.hh
  MashCorrectionPolicy.hh
  MeshScaledMassDensityPolicy.hh
  NonSymmetricSpecificThermalEnergyPolicy.hh
  PositionPolicy.hh
  PressurePolicy.hh
  SecondMomentHourglassControl.hh
  SoundSpeedPolicy.hh
  SpecificFromTotalThermalEnergyPolicy.hh
  SpecificThermalEnergyPolicy.hh
  SumVoronoiMassDensityPolicy.hh
  ThirdMomentHourglassControl.hh
  TotalHydro.hh
  TotalHydroInline.hh
  VolumePolicy.hh
  VoronoiHourglassControl.hh
  VoronoiMassDensityPolicy.hh
 )

set(Hydro_includes
  ${SPHERAL_SOURCE_DIR}
  ${SPHERAL_SOURCE_DIR}/Hydro
  ${EIGEN_LIB_LOCATION}
  )

blt_add_library(NAME        Spheral_Hydro
                HEADERS     ${Hydro_headers}
                SOURCES     ${Hydro_sources}
                DEPENDS_ON  ${spheral_depends}
                INCLUDES    ${Hydro_includes}
                )

install(TARGETS             Spheral_Hydro
        EXPORT              ${PROJECT_NAME}-targets
        DESTINATION         lib
        INCLUDES DESTINATION include
        )

install(FILES       ${Hydro_headers}
        DESTINATION include/Hydro
        ) 

list(APPEND SPHERAL_PYTHON_DEPENDS Spheral_Hydro)
set(SPHERAL_PYTHON_DEPENDS "${SPHERAL_PYTHON_DEPENDS}" PARENT_SCOPE)
