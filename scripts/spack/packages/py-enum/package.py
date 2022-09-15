# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class PyEnum(PythonPackage):
    """"""

    homepage = "https://pypi.org/project/enum/"
    pypi = "enum/enum-0.4.7.tar.gz" 

    #maintainers = ['mdavis36','jmikeowen']

    version('0.4.7', sha256='8c7cf3587eda51008bcc1eed99ea2c331ccd265c231dbaa95ec5258d3dc03100')

    extends('python@2.7:2.8', type=['build', 'run'])
