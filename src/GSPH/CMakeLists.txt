include_directories(.)

set(GSPH_inst
    GSPHHydroBase
    GSPHSpecificThermalEnergyPolicy
    Limiters/SlopeLimiterBase
    Limiters/VanLeerLimiter
    Limiters/SuperbeeLimiter
    Limiters/MinModLimiter
    Limiters/VanAlbaLimiter
    Limiters/OspreLimiter
    WaveSpeeds/WaveSpeedBase
    WaveSpeeds/AcousticWaveSpeed
    WaveSpeeds/DavisWaveSpeed
    WaveSpeeds/EinfeldtWaveSpeed
    RiemannSolvers/RiemannSolverBase)

set(GSPH_sources 
    GSPHFieldNames.cc)

set(GSPH_headers
    GSPHFieldNames.hh
    GSPHHydroBase.hh
    GSPHSpecificThermalEnergyPolicy.hh
    Limiters/SlopeLimiterBase.hh
    Limiters/VanLeerLimiter.hh
    Limiters/SuperbeeLimiter.hh
    Limiters/MinModLimiter.hh
    Limiters/VanAlbaLimiter.hh
    Limiters/OspreLimiter.hh
    WaveSpeeds/WaveSpeedBase.hh
    WaveSpeeds/AcousticWaveSpeed.hh
    WaveSpeeds/DavisWaveSpeed.hh
    WaveSpeeds/EinfeldtWaveSpeed.hh
    RiemannSolvers/RiemannSolverBase.hh)

add_subdirectory(Limiters)
add_subdirectory(WaveSpeeds)
add_subdirectory(RiemannSolvers)

instantiate(GSPH_inst GSPH_sources)

spheral_install_python_files( GSPHHydros.py)

spheral_add_cxx_library(GSPH)

