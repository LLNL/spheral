# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyGnuplot(PythonPackage):
    """Gnuplot.py is a Python package that allows you to create graphs from
       within Python using the gnuplot plotting program."""
    homepage = "https://pypi.org/project/py-gnuplot"
    url      = "https://files.pythonhosted.org/packages/1e/3f/2da7ee9232f8102ccbdef80b681d98ab286edd11b2320632f943c7828899/py-gnuplot-1.1.8.tar.gz"
    version('1.1.8', sha256='9c1404de6c27c76a5f43418a04c76c7706eb5238ba89781babae14280bfd1ada')

    # pip silently replaces distutils with setuptools
    depends_on('py-setuptools', type='build')
    depends_on('py-numpy', type=('build', 'run'))
