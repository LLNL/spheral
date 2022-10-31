# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyPipreqs(PythonPackage):
    """"""

    homepage = "https://pypi.org/project/pipreqs/"
    pypi = "pipreqs/pipreqs-0.4.10.tar.gz" 

    #maintainers = ['mdavis36','jmikeowen']

    version('0.4.10', sha256='9e351d644b28b98d7386b046a73806cbb3bb66b23a30e74feeb95ed9571db939')

    extends('python@2.7:2.8', type=['build', 'run'])
