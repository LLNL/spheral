include_directories(.)
set(SmoothingScale_inst
  SmoothingScaleBase
  SPHSmoothingScale
  ASPHSmoothingScale
  IncrementASPHHtensor
)

set(SmoothingScale_sources
  polySecondMoment.cc
)

instantiate(SmoothingScale_inst SmoothingScale_sources)

set(SmoothingScale_headers
  SmoothingScaleBase.hh
  SmoothingScaleBaseInline.hh
  FixedSmoothingScale.hh
  SPHSmoothingScale.hh
  ASPHSmoothingScale.hh
  ASPHSmoothingScaleUserFilter.hh
  polySecondMoment.hh
)

spheral_add_obj_library(SmoothingScale SPHERAL_OBJ_LIBS)
