########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             #
# Copyright (C) 2009 Sergio Galindo                                    #
# Copyright (C) 2013 William Oquendo                                   #
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
CMAKE_MINIMUM_REQUIRED (VERSION 3.28)

# --- Force the correct CUDA compiler and architecture ---
set(CMAKE_CUDA_COMPILER "/usr/local/cuda-13.0/bin/nvcc" CACHE FILEPATH "CUDA compiler" FORCE)
set(CMAKE_CUDA_FLAGS "-arch=sm_86" CACHE STRING "CUDA flags" FORCE)
set(CMAKE_CUDA_ARCHITECTURES "86" CACHE STRING "CUDA architectures" FORCE)

PROJECT (MECHSYSSANDBOX)

INCLUDE(${MECHSYSSANDBOX_SOURCE_DIR}/UseMECHSYS.cmake)


SET(PROGS
   tmpm01
   tmpm02
   tmpm03
)

SET(CUDAPROGS
    tmpm_cu_01
    tmpm_cu_02
    tmpm_cu_03
)

FOREACH(var ${PROGS})
    ADD_EXECUTABLE        (${var} "${var}.cpp")
    TARGET_LINK_LIBRARIES (${var} ${LIBS})
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")
ENDFOREACH(var)

FOREACH(var ${CUDAPROGS})
    ADD_EXECUTABLE        (${var} "${var}.cu")
    TARGET_LINK_LIBRARIES (${var} ${CUDA_LIBRARIES})
    # Remove any explicit -arch flag; rely on global CMAKE_CUDA_FLAGS
    SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS} -DUSE_CUDA" LINK_FLAGS "${LFLAGS}")
    # No need for target CUDA_ARCHITECTURES; global setting handles it
ENDFOREACH(var)
