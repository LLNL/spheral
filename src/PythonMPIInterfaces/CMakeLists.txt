if (ENABLE_MPI)
  set(PythonMPIInterfaces_python_targets
    mpi_mpi4py.py
    )
else()
  set(PythonMPIInterfaces_python_targets
    fakempi.py
    )
endif()

install(FILES       ${PythonMPIInterfaces_python_targets}
        DESTINATION Spheral
        RENAME      mpi.py
        )
