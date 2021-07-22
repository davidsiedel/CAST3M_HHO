# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# Maintainer comments:
# 18/12/2018: fix python detection

from spack import *


class Hho(CMakePackage):
    """
    HHO library for Cast3M
    """

    # homepage = "http://tfel.sourceforge.net"
    # url      = "https://github.com/thelfer/tfel/archive/TFEL-3.4.0.tar.gz"
    git      = "https://github.com/davidsiedel/CAST3M_HHO.git"
    maintainers = ['davidsiedel']

    # development branches
    version("main", branch="main")

    depends_on('eigen@3.3.8', type=('build', 'run'))

    extends('python', when='+python_bindings')

    def cmake_args(self):

        args = []

        args.append("-DCMAKE_BUILD_TYPE=Release")
        args.append("-DCMAKE_FIND_DEBUG_MODE=ON")
        args.append("-DTEST_MODE=test")

        return args