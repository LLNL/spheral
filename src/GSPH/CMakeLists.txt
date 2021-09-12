include_directories(.)

set(GSPH_inst
    GSPHHydroBase
    GSPHSpecificThermalEnergyPolicy
    Limiters/LimiterBase
    Limiters/VanLeerLimiter
    Limiters/SuperbeeLimiter
    Limiters/MinModLimiter
    Limiters/VanAlbaLimiter
    Limiters/OspreLimiter
    Limiters/PearlLimiter
    Limiters/BarthJespersenMinLimiter
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
    GSPHFieldNames.hh
    GSPHHydroBase.hh
    GSPHSpecificThermalEnergyPolicy.hh
    Limiters/LimiterBase.hh
    Limiters/VanLeerLimiter.hh
    Limiters/SuperbeeLimiter.hh
    Limiters/MinModLimiter.hh
    Limiters/VanAlbaLimiter.hh
    Limiters/OspreLimiter.hh
    Limiters/PearlLimiter.hh
    Limiters/BarthJespersenMinLimiter.hh
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

instantiate(GSPH_inst GSPH_sources)

spheral_install_python_files( GSPHHydros.py)

spheral_add_cxx_library(GSPH)

