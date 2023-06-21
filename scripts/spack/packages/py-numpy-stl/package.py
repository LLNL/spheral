# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyNumpyStl(PythonPackage):
    """"""

    homepage = "https://pypi.org/project/numpy-stl/"
    pypi = "numpy-stl/numpy-stl-2.11.2.tar.gz" 

    #maintainers = ['mdavis36','jmikeowen']

    version('2.11.2', sha256='192556df794b9ef0c1333fd5f034e4a3905d63f52345a0cc1e359045670e34b6')

    extends('python@2.7:2.8', type=['build', 'run'])
    depends_on('py-enum34', type='build')
    depends_on('py-python-utils@2.4.0', type='build')
    depends_on("py-setuptools", type="build")
