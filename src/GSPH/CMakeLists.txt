include_directories(.)

set(GSPH_sources 
    GSPHFieldNames.cc)

# Add explicit instantiations
list(APPEND GSPH_sources
     computeSPHVolumeInst.cc
     computeSumVolumeInst.cc
     computeMFMDensityInst.cc
     initializeGradientsInst.cc
     GenericRiemannHydroInst.cc
     GSPHHydroBaseInst.cc
     MFMHydroBaseInst.cc
     MFVHydroBaseInst.cc
     Policies/MassFluxPolicyInst.cc
     Policies/ReplaceWithRatioPolicyInst.cc
     Policies/MFVIncrementVelocityPolicyInst.cc
     Policies/MFVIncrementSpecificThermalEnergyPolicyInst.cc
     Policies/CompatibleMFVSpecificThermalEnergyPolicyInst.cc
     Limiters/LimiterBaseInst.cc
     Limiters/VanLeerLimiterInst.cc
     Limiters/SuperbeeLimiterInst.cc
     Limiters/MinModLimiterInst.cc
     Limiters/VanAlbaLimiterInst.cc
     Limiters/OspreLimiterInst.cc
     Limiters/BarthJespersenLimiterInst.cc
     WaveSpeeds/WaveSpeedBaseInst.cc
     WaveSpeeds/AcousticWaveSpeedInst.cc
     WaveSpeeds/DavisWaveSpeedInst.cc
     WaveSpeeds/EinfeldtWaveSpeedInst.cc
     RiemannSolvers/RiemannSolverBaseInst.cc
     RiemannSolvers/HLLCInst.cc
     RiemannSolvers/SecondOrderArtificialViscosityInst.cc
     )

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
    Policies/ReplaceWithRatioPolicy.hh
    Policies/MFVIncrementVelocityPolicy.hh
    Policies/MFVIncrementSpecificThermalEnergyPolicy.hh
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

spheral_install_python_files(GSPHHydros.py)

spheral_add_obj_library(GSPH SPHERAL_OBJ_LIBS)

