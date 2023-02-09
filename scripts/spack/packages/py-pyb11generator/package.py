# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyPyb11generator(PythonPackage):
    """PYB11Generator is a python based code generator that creates pybind11 code for binding C++ libraries as extensions in Python."""

    homepage = "https://pypi.org/project/PYB11Generator/"
    git      = "https://github.com/LLNL/PYB11Generator.git"
    #pypi     = "PYB11Generator/PYB11Generator-2.0.2.tar.gz" 

    maintainers = ['mdavis36','jmikeowen']

    version('2.0.2', tag='2.0.2', submodules=True)
    #version('2.0.2', sha256='671941aaaf872a5837ebd73bc9c6cd43f4fdba39225cd153c128411478b7aa17')

    extends('python@3:', type=['build', 'run'])
    #depends_on('py-pybind11', type=('build', 'run'))   # Getting from PYB11Generator submodule
    depends_on('py-setuptools', type=('build', 'run'))
    depends_on('py-decorator', type=('build', 'run'))
