include_directories(.)
set(SmoothingScale_inst
    SmoothingScaleBase
    SPHSmoothingScale
    ASPHSmoothingScale
   )

set(SmoothingScale_sources )

instantiate(SmoothingScale_inst SmoothingScale_sources)

set(SmoothingScale_headers
    SmoothingScaleBase.hh
    SmoothingScaleBaseInline.hh
    FixedSmoothingScale.hh
    SPHSmoothingScale.hh
    ASPHSmoothingScale.hh
    ASPHSmoothingScaleUserFilter.hh
    )

spheral_add_obj_library(SmoothingScale SPHERAL_OBJ_LIBS)
