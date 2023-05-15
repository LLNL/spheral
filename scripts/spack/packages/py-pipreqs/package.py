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

    version('0.4.11', sha256='c793b4e147ac437871b3a962c5ce467e129c859ece5ba79aca83c20f4d9c3aef')

    extends('python@3:', type=['build', 'run'])
    depends_on("py-setuptools", type="build")
