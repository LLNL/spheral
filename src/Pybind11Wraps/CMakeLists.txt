set(Pybind11Wraps_python_targets
  SpheralPickle.py
  PolyhedronFileUtilities.py
  )

install(FILES       ${Pybind11Wraps_python_targets}
        DESTINATION .
        )
