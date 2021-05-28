include_directories(.)
set(DeviceTestLib2_inst
    )


set(DeviceTestLib2_sources 
    DeviceTest2.cc)

instantiate(DeviceTestLib2_inst DeviceTestLib2_sources)

set(DeviceTestLib2_headers
    DeviceTest2.hh)

spheral_install_python_files(
    )

spheral_add_cxx_library(DeviceTestLib2)

