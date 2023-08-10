include_directories(.)

set(GSPH_inst
    computeSPHVolume
    computeSumVolume
    computeMFMDensity
    initializeGradients
    GenericRiemannHydro
    GSPHHydroBase
    MFMHydroBase
    MFVHydroBase
    Policies/MassFluxPolicy
    Policies/ALEPositionPolicy
    Policies/PureReplaceWithStateFieldList
    Policies/PureReplaceFieldList
    Policies/ReplaceWithRatioPolicy
    Policies/IncrementSpecificFromTotalPolicy
    Policies/CompatibleMFVSpecificThermalEnergyPolicy
    Limiters/LimiterBase
    Limiters/VanLeerLimiter
    Limiters/SuperbeeLimiter
    Limiters/MinModLimiter
    Limiters/VanAlbaLimiter
    Limiters/OspreLimiter
    Limiters/BarthJespersenLimiter
    WaveSpeeds/WaveSpeedBase
    WaveSpeeds/AcousticWaveSpeed
    WaveSpeeds/DavisWaveSpeed
    WaveSpeeds/EinfeldtWaveSpeed
    RiemannSolvers/RiemannSolverBase
    RiemannSolvers/HLLC
    RiemannSolvers/SecondOrderArtificialViscosity)

set(GSPH_sources 
    GSPHFieldNames.cc)

set(GSPH_headers
    computeSPHVolume.hh
    computeSumVolume.hh
    computeMFMDensity.hh
    initializeGradients.hh
    GSPHFieldNames.hh
    GenericRiemannHydro.hh
    GSPHHydroBase.hh
    MFMHydroBase.hh
    MFVHydroBase.hh
    Policies/MassFluxPolicy.hh
    Policies/ALEPositionPolicy.hh
    Policies/PureReplaceWithStateFieldList.hh
    Policies/PureReplaceFieldList.hh
    Policies/ReplaceWithRatioPolicy.hh
    Policies/IncrementSpecificFromTotalPolicy.hh
    Policies/CompatibleMFVSpecificThermalEnergyPolicy.hh
    Limiters/LimiterBase.hh
    Limiters/VanLeerLimiter.hh
    Limiters/SuperbeeLimiter.hh
    Limiters/MinModLimiter.hh
    Limiters/VanAlbaLimiter.hh
    Limiters/OspreLimiter.hh
    Limiters/BarthJespersenLimiter.hh
    WaveSpeeds/WaveSpeedBase.hh
    WaveSpeeds/AcousticWaveSpeed.hh
    WaveSpeeds/DavisWaveSpeed.hh
    WaveSpeeds/EinfeldtWaveSpeed.hh
    RiemannSolvers/RiemannSolverBase.hh
    RiemannSolvers/HLLC.hh
    RiemannSolvers/SecondOrderArtificialViscosity.hh)

add_subdirectory(Limiters)
add_subdirectory(WaveSpeeds)
add_subdirectory(RiemannSolvers)
add_subdirectory(Policies)

instantiate(GSPH_inst GSPH_sources)

spheral_install_python_files(GSPHHydros.py)

spheral_add_obj_library(GSPH)

