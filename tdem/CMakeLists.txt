########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2009 Sergio Galindo                                    #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

SET(EXES
  nonconvex
  ttt
  test_periodic
  test_distances
  test_tetrahedra
  test_domain
  test_dynamics
  test_mesh
  test_read
  test_write
  test_beam
  cylinder_test
  brazil_test
  test_01
  test_02
  GSD
)

SET(CUDAPROGS
    test_cu_01
    test_cu_02
    test_cu_03
    test_cu_04
)

FOREACH(var ${EXES})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

FOREACH(var ${CUDAPROGS})
    ADD_EXECUTABLE        (${var} "${var}.cu")
    TARGET_LINK_LIBRARIES (${var} ${CUDA_LIBRARIES})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS} -DUSE_CUDA" LINK_FLAGS "${LFLAGS}")
    #SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS} -DUSE_CUDA -arch=sm_60" LINK_FLAGS "${LFLAGS}")
    SET_TARGET_PROPERTIES (${var} PROPERTIES CUDA_ARCHITECTURES "native")
ENDFOREACH(var)
