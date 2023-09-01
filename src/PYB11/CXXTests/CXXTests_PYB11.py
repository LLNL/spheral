"""
Spheral CXXTests module.

This module provide thin front-end wrappers for C++ in Spheral.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#PYB11includes = ['"CXXTests/DeviceTestLib/DeviceTest.hh"']

PYB11namespaces = ["Spheral"]

def launchCaller(a = "int",
                 b = "int"):
  return "int"

