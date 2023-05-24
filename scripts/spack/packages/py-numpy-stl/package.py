# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyNumpyStl(PythonPackage):
    """"""

    homepage = "https://pypi.org/project/numpy-stl/"
    pypi = "numpy-stl/numpy-stl-3.0.0.tar.gz" 

    version('3.0.0', sha256='578b78eacb0529ac9aba2f17dcc363d58c7c3c5708710c18f8c1e9965f2e81ac')

    extends('python@3:', type=['build', 'run'])
    depends_on('py-setuptools', type='build')
