# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyPyb11generator(PythonPackage):
    """PYB11Generator is a python based code generator that creates pybind11 code for binding C++ libraries as extensions in Python."""

    homepage = "https://pypi.org/project/PYB11Generator/"
    pypi = "PYB11Generator/PYB11Generator-1.0.12.tar.gz" 

    maintainers = ['mdavis36','jmikeowen']

    version('1.0.12', sha256='0a0988e705aebf050180170b25b57e4a6c7652d2099f2ef180eed63c4712d91c')

    extends('python@2.7:2.8', type=['build', 'run'])
    depends_on('py-setuptools', type=('build', 'run'))
    depends_on('py-decorator', type=('build', 'run'))
