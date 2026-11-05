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
  test_2spheres
  PID_v2
  TX_PID_periodic
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
  test_saveload
  GSD
)

SET(CUDAPROGS
    test_2spheres_gpu
    test_cu_01
    test_cu_02
    test_cu_03
    test_cu_04
)

# GPU build of the periodic PID driver: same base name as the CPU driver
# (TX_PID_periodic.cu), so the executable is suffixed with _gpu.
ADD_EXECUTABLE        (TX_PID_periodic_gpu "TX_PID_periodic.cu")
TARGET_LINK_LIBRARIES (TX_PID_periodic_gpu ${CUDA_LIBRARIES})
SET_TARGET_PROPERTIES (TX_PID_periodic_gpu PROPERTIES COMPILE_FLAGS "${FLAGS} -DUSE_CUDA" LINK_FLAGS "${LFLAGS}")
IF(NOT MECHSYS_CUDA_ARCH_VIA_FLAGS)
    SET_TARGET_PROPERTIES (TX_PID_periodic_gpu PROPERTIES CUDA_ARCHITECTURES "${MECHSYS_CUDA_ARCH}")
ENDIF()

# Constant-volume Lees-Edwards shearbox: a production-style example driven
# entirely through Domain's function pointers (pForceVV, pReset, pExtraParams)
# without needing any Lees-Edwards support inside the library itself.
LIST(APPEND CUDAPROGS shearbox_LE_const_vol)

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
    # MECHSYS_CUDA_ARCH_VIA_FLAGS is set when the CMake in use is too old for
    # this property; the root CMakeLists.txt already passed -gencode directly.
    IF(NOT MECHSYS_CUDA_ARCH_VIA_FLAGS)
        SET_TARGET_PROPERTIES (${var} PROPERTIES CUDA_ARCHITECTURES "${MECHSYS_CUDA_ARCH}")
    ENDIF()
ENDFOREACH(var)
