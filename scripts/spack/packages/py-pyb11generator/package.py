# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install pyb11generator
#
# You can edit this file again by typing:
#
#     spack edit pyb11generator
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class PyPyb11generator(PythonPackage):
    """FIXME: Put a proper description of your package here."""

    homepage = "https://pypi.org/project/PYB11Generator/"
    pypi = "PYB11Generator/PYB11Generator-1.0.12.tar.gz" 

    maintainers = ['mdavis36','jmikeowen']

    version('1.0.12', sha256='0a0988e705aebf050180170b25b57e4a6c7652d2099f2ef180eed63c4712d91c')

    extends('python@2.7:2.8', type=['build', 'run'])
    depends_on('py-setuptools', type=('build', 'run'))
    depends_on('py-decorator')
