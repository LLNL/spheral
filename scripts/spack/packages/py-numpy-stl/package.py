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
#     spack install py-numpy-stl
#
# You can edit this file again by typing:
#
#     spack edit py-numpy-stl
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack import *


class PyNumpyStl(PythonPackage):
    """FIXME: Put a proper description of your package here."""

    # FIXME: Add a proper url for your package's homepage here.
    homepage      = "https://pypi.org/project/numpy-stl/"
    pypi = "numpy-stl/numpy-stl-2.11.2.tar.gz"

    version('2.11.2', 'd625e8c11a6cfb475d9c33781593e830')

    extends('python@2.7:2.8', type=['build','run'])
    extends('py-setuptools', type=['build','run'])
    # FIXME: Add dependencies if required.
    # depends_on('foo')
