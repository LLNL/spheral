include_directories(.)

set(GSPH_inst
    computeSPHVolume
    computeSumVolume
    computeMFMDensity
    GenericRiemannHydro
    GSPHHydroBase
    MFMHydroBase
    Policies/PureReplaceFieldList
    Policies/ReplaceWithRatioPolicy
    Limiters/LimiterBase
    Limiters/VanLeerLimiter
    Limiters/SuperbeeLimiter
    Limiters/MinModLimiter
    Limiters/VanAlbaLimiter
    Limiters/OspreLimiter
    WaveSpeeds/WaveSpeedBase
    WaveSpeeds/AcousticWaveSpeed
    WaveSpeeds/DavisWaveSpeed
    WaveSpeeds/EinfeldtWaveSpeed
    RiemannSolvers/RiemannSolverBase
    RiemannSolvers/HLLC
    RiemannSolvers/GHLLC)

set(GSPH_sources 
    GSPHFieldNames.cc)

set(GSPH_headers
    computeSPHVolume.hh
    computeSumVolume.hh
    computeMFMDensity.hh
    GSPHFieldNames.hh
    GenericRiemannHydro.hh
    GSPHHydroBase.hh
    MFMHydroBase.hh
    Policies/PureReplaceFieldList.hh
    Policies/ReplaceWithRatioPolicy.hh
    Limiters/LimiterBase.hh
    Limiters/VanLeerLimiter.hh
    Limiters/SuperbeeLimiter.hh
    Limiters/MinModLimiter.hh
    Limiters/VanAlbaLimiter.hh
    Limiters/OspreLimiter.hh
    WaveSpeeds/WaveSpeedBase.hh
    WaveSpeeds/AcousticWaveSpeed.hh
    WaveSpeeds/DavisWaveSpeed.hh
    WaveSpeeds/EinfeldtWaveSpeed.hh
    RiemannSolvers/RiemannSolverBase.hh
    RiemannSolvers/HLLC.hh
    RiemannSolvers/GHLLC.hh)

add_subdirectory(Limiters)
add_subdirectory(WaveSpeeds)
add_subdirectory(RiemannSolvers)
add_subdirectory(Policies)

instantiate(GSPH_inst GSPH_sources)

spheral_install_python_files(GSPHHydros.py)

spheral_add_cxx_library(GSPH)

