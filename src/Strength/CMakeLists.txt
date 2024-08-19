include_directories(.)

set(Strength_sources
    SolidFieldNames.cc
    )

# Add explicit instantiations
list(APPEND Strength_sources
     BulkModulusPolicyInst.cc
     DeviatoricStressPolicyInst.cc
     MeltEnergyPolicyInst.cc
     PlasticStrainPolicyInst.cc
     ShearModulusPolicyInst.cc
     YieldStrengthPolicyInst.cc
     )

set(Strength_headers
    BulkModulusPolicy.hh
    DeviatoricStressPolicy.hh
    MeltEnergyPolicy.hh
    PlasticStrainPolicy.hh
    ShearModulusPolicy.hh
    SolidFieldNames.hh
    YieldStrengthPolicy.hh
    effectiveKernelVolume.hh
    )

spheral_add_obj_library(Strength SPHERAL_OBJ_LIBS)

