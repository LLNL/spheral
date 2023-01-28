# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyPyb11generator(PythonPackage):
    """PYB11Generator is a python based code generator that creates pybind11 code for binding C++ libraries as extensions in Python."""

    homepage = "https://pypi.org/project/PYB11Generator/"
    pypi = "PYB11Generator/PYB11Generator-2.0.1.tar.gz" 

    maintainers = ['mdavis36','jmikeowen']

    version('2.0.1', sha256='a283ffccb2a4a0cb0bc4ae6470bc710cb60869edfbff75fee3631401d3b35dc9')

    extends('python@3:', type=['build', 'run'])
    depends_on('py-pybind11', type=('build', 'run'))
    depends_on('py-setuptools', type=('build', 'run'))
    depends_on('py-decorator', type=('build', 'run'))
